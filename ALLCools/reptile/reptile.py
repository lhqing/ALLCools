import pathlib
import numpy as np
import pandas as pd
import xarray as xr
import pybedtools
import dask
import subprocess
from collections import defaultdict
from ALLCools.mcds import RegionDS
import pyBigWig
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.model_selection import train_test_split
from ALLCools.mcds.utilities import write_ordered_chunks, update_dataset_config
import joblib



def _create_train_region_ds(reptile):
    # train regions
    train_regions_bed = pybedtools.BedTool(reptile.train_regions).sort(
        g=reptile.chrom_size_path
    )

    # train region labels
    train_label = pd.read_csv(
        reptile.train_region_labels, sep="\t", index_col=0, squeeze=True
    )
    train_label.index.name = "train-region"
    train_regions_bed_df = train_regions_bed.to_dataframe()

    # train RegionDS
    train_region_ds = RegionDS.from_bed(
        train_regions_bed_df,
        chrom_size_path=reptile.chrom_size_path,
        location=reptile.output_path,
        region_dim="train-region",
    )

    train_region_ds.coords["train-region_label"] = train_label
    train_region_ds.save()
    return train_regions_bed, train_label


def _create_train_dmr_ds(reptile, train_regions_bed, train_label):
    # total DMRs
    dmr_regions_bed = pybedtools.BedTool(reptile.dmr_regions).sort(
        g=reptile.chrom_size_path
    )
    dmr_regions_bed_df = dmr_regions_bed.to_dataframe()

    # train DMRs and train DMR labels
    train_dmr = train_regions_bed.map(dmr_regions_bed, c=4, o="collapse").to_dataframe()

    dmr_label = defaultdict(list)
    for _, row in train_dmr.iterrows():
        *_, train_region, dmrs = row
        if dmrs == ".":
            continue
        dmrs = dmrs.split(",")
        for dmr in dmrs:
            dmr_label[dmr].append(train_label[train_region])

    # some DMR might have multiple labels
    consistent_dmr_label = {}
    for dmr, dmr_labels in dmr_label.items():
        if (len(dmr_labels) == 1) or (len(set(dmr_labels)) == 1):
            consistent_dmr_label[dmr] = dmr_labels[0]
        else:
            # dmr has in consistent label
            continue
    dmr_label = pd.Series(consistent_dmr_label)
    dmr_label.index.name = "train-dmr"

    train_dmr_regions_bed_df = (
        dmr_regions_bed_df.set_index("name")
            .loc[dmr_label.index]
            .reset_index()
            .iloc[:, [1, 2, 3, 0]]
    )

    # train DMR RegionDS
    train_dmr_ds = RegionDS.from_bed(
        train_dmr_regions_bed_df,
        chrom_size_path=reptile.chrom_size_path,
        location=reptile.output_path,
        region_dim="train-dmr",
    )

    train_dmr_ds.coords["train-dmr_label"] = dmr_label
    train_dmr_ds.save()
    return dmr_regions_bed_df


def _create_query_region_ds(reptile):
    pybedtools.BedTool().makewindows(
        g=reptile.chrom_size_path, s=reptile.step_size, w=reptile.window_size
    ).saveas(f"{reptile.output_path}/query_region.bed")

    query_region_ds = RegionDS.from_bed(
        f"{reptile.output_path}/query_region.bed",
        chrom_size_path=reptile.chrom_size_path,
        location=reptile.output_path,
        region_dim="query-region",
    )

    subprocess.run(f"rm -f {reptile.output_path}/query_region.bed", shell=True)
    query_region_ds.save()
    return


def _create_query_dmr_ds(reptile, dmr_regions_bed_df):
    query_dmr_ds = RegionDS.from_bed(
        dmr_regions_bed_df,
        chrom_size_path=reptile.chrom_size_path,
        location=reptile.output_path,
        region_dim="query-dmr",
    )
    query_dmr_ds.save()
    return


def _get_data_and_label(region_ds, modalities, sample, fillna_by_zero_list):
    data = {}
    for modality in modalities:
        da = region_ds[f"{region_ds.region_dim}_{modality}_da"]
        if modality in fillna_by_zero_list:
            da = da.fillna(0)

        sample_value = da.sel({modality: sample}).to_pandas()

        modality_mean = da.coords[f"{modality}_mean"]
        modality_std = da.coords[f"{modality}_std"]
        sample_dev = sample_value - modality_mean

        features = {
            f"{modality}_value": sample_value,
            f"{modality}_dev": sample_dev,
            f"{modality}_overall_std": modality_std,
        }
        data.update(features)

    data = pd.DataFrame(data)
    if region_ds.region_dim.startswith("train"):
        label = region_ds[f"{region_ds.region_dim}_label"].to_pandas()
    else:
        label = None
    return data, label


def _predict_sample(
        region_ds_path,
        region_dim,
        modalities,
        fillna_by_zero,
        sample,
        output_path,
        mask_cutoff=0.3,
        chunk_size=100000,
):
    # set dask scheduler to allow multiprocessing
    with dask.config.set(scheduler="sync"):
        if region_dim == "query-region":
            model = joblib.load(f"{region_ds_path}/model/train-region_model.lib")
        elif region_dim == "query-dmr":
            model = joblib.load(f"{region_ds_path}/model/train-dmr_model.lib")
        else:
            raise ValueError(
                f'Only accept ["query-region", "query-dmr"], got {region_dim}'
            )

        # for query, we don't filter nan, but drop the nan value in final data table
        region_ds = RegionDS.open(region_ds_path, region_dim=region_dim)
        region_ids = region_ds.get_index(region_dim)

        total_proba = []
        for chunk_start in range(0, region_ids.size, chunk_size):
            use_regions = region_ids[chunk_start: chunk_start + chunk_size]
            _region_ds = region_ds.sel({region_dim: use_regions})
            data, _ = _get_data_and_label(
                region_ds=_region_ds,
                modalities=modalities,
                sample=sample,
                fillna_by_zero_list=fillna_by_zero,
            )
            # before dropna, save the index
            total_index = data.index.copy()

            # sample specific NaN drop
            data.dropna(inplace=True)

            # predict
            proba = model.predict_proba(data.astype(np.float64))
            enhancer_proba = pd.Series(proba[:, 1], index=data.index).reindex(
                total_index
            )
            # NA value has 0 proba
            enhancer_proba.fillna(0, inplace=True)
            total_proba.append(enhancer_proba)

        total_proba = pd.DataFrame({sample: pd.concat(total_proba).astype(np.float16)})
        # mask small values
        total_proba[total_proba < mask_cutoff] = 0
        total_proba.index.name = region_ds.region_dim
        total_proba.columns.name = "sample"

        total_proba = xr.Dataset({f"{region_dim}_prediction": total_proba})
        RegionDS(total_proba).to_zarr(output_path, mode="w")
        return output_path


class REPTILE:
    def __init__(
            self,
            output_path,
            train_regions,
            dmr_regions,
            train_region_labels,
            train_sample,
            bigwig_table,
            chrom_size_path,
            window_size=2000,
            step_size=200,
            dmr_slop=150,
            fillna_by_zero=None,
    ):
        try:
            from tpot import TPOTClassifier
        except ImportError:
            raise ImportError('Got tpot import error, install the tpot package: \n'
                              'conda install -c conda-forge tpot xgboost dask dask-ml scikit-mdr skrebate \n'
                              'See also https://epistasislab.github.io/tpot/')
        self.output_path = output_path
        pathlib.Path(self.output_path).mkdir(exist_ok=True)
        # for reptile model
        self.model_dir = f"{self.output_path}/model"
        pathlib.Path(self.model_dir).mkdir(exist_ok=True)
        # for final prediction results
        self.bigwig_dir = f"{self.output_path}/bigwig"
        pathlib.Path(self.bigwig_dir).mkdir(exist_ok=True)

        self.train_regions = train_regions
        self.dmr_regions = dmr_regions
        self.train_region_labels = train_region_labels
        self.train_sample = train_sample
        self.chrom_size_path = chrom_size_path
        self.window_size = window_size
        self.step_size = step_size
        self.dmr_slop = dmr_slop

        if fillna_by_zero is None:
            fillna_by_zero = []
        self.fillna_by_zero = fillna_by_zero

        self.bigwig_table: pd.DataFrame
        if isinstance(bigwig_table, (str, pathlib.PosixPath)):
            bigwig_table = str(bigwig_table)
            if bigwig_table.endswith("tsv"):
                sep = "\t"
            else:
                sep = ","
            self.bigwig_table = pd.read_csv(bigwig_table, sep=sep, index_col=0)
        else:
            self.bigwig_table = bigwig_table

        self.modalities = self.bigwig_table.columns
        print(f"Got {self.modalities.size} modalities from bigwig_table: ", end="")
        print(", ".join(self.modalities))
        if train_sample not in self.bigwig_table.index:
            raise KeyError(
                f"train_sample {train_sample} does not exist in the index of bigwig_table"
            )
        self.samples = pd.Index(
            [s for s in self.bigwig_table.index if s != train_sample]
        )
        print("Training sample:", self.train_sample)
        print("Other samples:", ", ".join(self.samples))

        # four RegionDS
        self._train_region_ds = None
        self._train_dmr_ds = None
        self._query_region_ds = None
        self._query_dmr_ds = None
        self.region_dims = ["train-region", "train-dmr", "query-region", "query-dmr"]

        # two models
        self._region_model = None
        self._dmr_model = None

        # prediction results
        self._region_prediction = None
        self._dmr_prediction = None

        # check if default da exist, otherwise generate
        try:
            assert pathlib.Path(f"{self.output_path}/train-region").exists()
            assert pathlib.Path(f"{self.output_path}/train-dmr").exists()
            assert pathlib.Path(f"{self.output_path}/query-region").exists()
            assert pathlib.Path(f"{self.output_path}/query-dmr").exists()
        except AssertionError:
            self.generate_region_ds()
        return

    def generate_region_ds(self):
        # RegionDS for training
        # step 1. Create training region RegionDS
        train_regions_bed, train_label = _create_train_region_ds(self)
        # step 2. Create DMR RegionDS overlapping with training regions
        dmr_regions_bed_df = _create_train_dmr_ds(
            self, train_regions_bed=train_regions_bed, train_label=train_label
        )

        # RegionDS for query
        # step 3. Create all DMR RegionDS for query
        _create_query_dmr_ds(self, dmr_regions_bed_df=dmr_regions_bed_df)
        # step 4. Create all genome window (self.window_size, self.step) RegionDS for query
        _create_query_region_ds(self)
        return

    @property
    def train_region_ds(self):
        if self._train_region_ds is None:
            self._train_region_ds = RegionDS.open(
                self.output_path, region_dim="train-region", engine="zarr"
            )
        return self._train_region_ds

    @property
    def train_dmr_ds(self):
        if self._train_dmr_ds is None:
            self._train_dmr_ds = RegionDS.open(
                self.output_path, region_dim="train-dmr", engine="zarr"
            )
        return self._train_dmr_ds

    @property
    def query_region_ds(self):
        if self._query_region_ds is None:
            self._query_region_ds = RegionDS.open(
                self.output_path, region_dim="query-region", engine="zarr"
            )
        return self._query_region_ds

    @property
    def query_dmr_ds(self):
        if self._query_dmr_ds is None:
            self._query_dmr_ds = RegionDS.open(self.output_path, region_dim="query-dmr")
        return self._query_dmr_ds

    @property
    def region_model(self):
        if self._region_model is None:
            try:
                self._region_model = joblib.load(
                    f"{self.model_dir}/train-region_model.lib"
                )
            except FileNotFoundError:
                print(
                    "Region model not found, make sure you train the model before prediction"
                )
        return self._region_model

    @property
    def dmr_model(self):
        if self._dmr_model is None:
            try:
                self._dmr_model = joblib.load(f"{self.model_dir}/train-dmr_model.lib")
            except FileNotFoundError:
                print(
                    "DNR model not found, make sure you train the model before prediction"
                )
        return self._dmr_model

    def _validate_region_name(self, name):
        if name not in self.region_dims:
            raise KeyError(f"Only accept {self.region_dims}, got {name}")

    def annotate_by_bigwigs(self, name, slop, cpu, redo=False):
        attribute_name = f"{name.replace('-', '_')}_ds"
        region_ds: RegionDS = self.__getattribute__(attribute_name)
        if region_ds.attrs.get("reptile_annotate") and not redo:
            print(f"{name} already annotated")
            return

        for modality, bigwig_paths in self.bigwig_table.items():
            print(f"Annotating regions with {modality} BigWigs.")
            region_ds.annotate_by_bigwigs(
                bigwig_paths.to_dict(),
                dim=modality,
                slop=slop,  # None for region, 150 for DMR
                value_type="mean",
                chunk_size="auto",
                dtype="float32",
                cpu=cpu,
            )
            da_name = f"{region_ds.region_dim}_{modality}_da"
            if modality in self.fillna_by_zero:
                region_ds[da_name] = region_ds[da_name].fillna(0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                region_ds.coords[f"{modality}_mean"] = (
                    region_ds[da_name]
                        .sel({modality: self.samples})
                        .mean(dim=modality, skipna=True)
                        .load()
                )
                region_ds.coords[f"{modality}_std"] = (
                    region_ds[da_name]
                        .sel({modality: self.samples})
                        .std(dim=modality, skipna=True)
                        .load()
                )
                region_ds.coords[f"{modality}_nan_count"] = (
                    np.isnan(region_ds[da_name].sel({modality: self.samples}))
                        .sum(dim=modality)
                        .load()
                )
        # save the annotated region_ds to self.{output_dir}/{region_dim}
        region_ds.attrs["reptile_annotate"] = True

        region_ds.save(mode="a")
        for modality in self.modalities:
            subprocess.run(f"rm -rf {self.output_path}/{name}_{modality}", shell=True)

        # need to re-open region ds otherwise the annotated da is missing
        # will reopen when access this attr
        self.__setattr__(f"_{attribute_name}", None)
        return

    def _filter_na_train(self, name, sample, max_na_rate=0.5):
        self._validate_region_name(name)
        attribute_name = f"{name.replace('-', '_')}_ds"
        region_ds: RegionDS = self.__getattribute__(attribute_name)
        for modality in self.modalities:
            if modality not in region_ds.dims:
                raise KeyError(
                    f"{name} RegionDS missing {modality} annotation, "
                    f"make sure you run REPTILE.annotate_by_bigwigs before training."
                )

        judges = []
        for data_var, da in region_ds.data_vars.items():
            modality = data_var.split("_")[1]
            if modality in self.fillna_by_zero:
                # for count based data, nan means 0 count, no need to filter
                continue

            # remove regions having no coverage in training sample
            nan_sample = np.isnan(da.sel({modality: sample})).to_pandas()
            # or having <25% samples covered in other samples
            nan_other = da.coords[f"{modality}_nan_count"].to_pandas()
            nan_other = nan_other > (len(self.samples) * max_na_rate)
            remove_feature = nan_sample | nan_other
            judges.append(remove_feature)

        if len(judges) > 0:
            remove_feature = np.any(judges, axis=0)
            region_ds = region_ds.sel({region_ds.region_dim: ~remove_feature})
            print(
                f"{(~remove_feature).sum()} features remained after filter NaN for {name}"
            )
        else:
            # this means all modality NaN filled by 0
            pass
        return region_ds

    def prepare_training_input(self, name):
        if name not in self.region_dims:
            raise ValueError(f'Only accept ["train-region", "train-dmr"], got {name}')

        # get the RegionDS that with NaN values filtered (see self._filter_na)
        region_ds = self._filter_na_train(
            name=name, sample=self.train_sample, max_na_rate=0.5
        )
        data, label = _get_data_and_label(
            region_ds=region_ds,
            modalities=self.modalities,
            sample=self.train_sample,
            fillna_by_zero_list=self.fillna_by_zero,
        )
        return data, label

    @staticmethod
    def auto_ml(
            data,
            label,
            output_path,
            train_size=0.75,
            random_state=42,
            cpu=1,
            tpot_generations=5,
            tpot_max_time_mins=60,
            **tpot_kwargs,
    ):
        print("Training model with these parameters:")
        print(f"cpu={cpu}")
        print(f"train_size={train_size}")
        print(f"random_state={random_state}")
        print(f"generations={tpot_generations}")
        print(f"max_time_mins={tpot_max_time_mins}")
        for k, v in tpot_kwargs.items():
            print(f"{k}={v}")

        x_train, x_test, y_train, y_test = train_test_split(
            data.astype(np.float64),
            label.astype(np.float64),
            train_size=train_size,
            test_size=1 - train_size,
            random_state=random_state,
        )

        from tpot import TPOTClassifier
        _tpot = TPOTClassifier(
            generations=tpot_generations,
            max_time_mins=tpot_max_time_mins,
            verbosity=2,
            n_jobs=cpu,
            random_state=random_state,
            **tpot_kwargs,
        )
        _tpot.fit(x_train, y_train)
        final_score = _tpot.score(x_test, y_test)
        print(f"Final hold-out testing data accuracy: {final_score:.4f}")
        print("Final pipeline:")
        print(_tpot.fitted_pipeline_)

        # save the model for prediction
        joblib.dump(_tpot.fitted_pipeline_, output_path)
        return

    def _train(self, region_dim, slop, cpu, **kwargs):
        # step 1: annotate RegionDS
        self.annotate_by_bigwigs(region_dim, slop=slop, cpu=cpu)

        # step 2: prepare data matrix and filter nan
        data, label = self.prepare_training_input(region_dim)

        # step 3: perform AutoML training and get the final model
        model_output_path = f"{self.model_dir}/{region_dim}_model.lib"
        # add cpu otherwise auto_ml don't have this parameter
        kwargs["cpu"] = cpu
        self.auto_ml(data=data, label=label, output_path=model_output_path, **kwargs)
        return model_output_path

    def train_region_model(self, slop=None, cpu=1, **kwargs):
        region_dim = "train-region"
        print("Training model for genomic regions.")
        model_path = self._train(region_dim, slop, cpu=cpu, **kwargs)
        self._region_model = joblib.load(model_path)
        return

    def train_dmr_model(self, slop=None, cpu=1, **kwargs):
        region_dim = "train-dmr"
        print("Training model for DMR regions.")
        if slop is None:
            slop = self.dmr_slop
            print(f"Using slop {slop} for DMRs")
        model_path = self._train(region_dim, slop=slop, cpu=cpu, **kwargs)
        self._dmr_model = joblib.load(model_path)
        return

    def fit(self, cpu=10, **kwargs):
        """Convenient function to train everything by default parameters"""
        self.train_region_model(cpu=cpu, **kwargs)
        self.train_dmr_model(cpu=cpu, **kwargs)
        return

    def _predict(self, region_dim, cpu, mask_cutoff):
        if region_dim == "query-region":
            model = self.region_model
            slop = None
        else:
            model = self.dmr_model
            slop = self.dmr_slop
        print(f"Perform {region_dim} prediction using fitted model:")
        print(model)

        # step 1: annotate RegionDS
        self.annotate_by_bigwigs(region_dim, slop=slop, cpu=cpu)

        # step 2: prepare input and run prediction
        chunk_dir = f"{self.output_path}/.{region_dim}_prediction_chunks"
        pathlib.Path(chunk_dir).mkdir(exist_ok=True)
        final_dir = f"{self.output_path}/{region_dim}_prediction"

        with ProcessPoolExecutor(cpu) as exe:
            futures = {}
            for chunk_i, sample in enumerate(self.samples):
                sample_output_path = f"{chunk_dir}/{sample}.zarr"
                future = exe.submit(
                    _predict_sample,
                    region_ds_path=self.output_path,
                    region_dim=region_dim,
                    modalities=self.modalities,
                    fillna_by_zero=self.fillna_by_zero,
                    sample=sample,
                    output_path=sample_output_path,
                    chunk_size=500000,
                    mask_cutoff=mask_cutoff,
                )
                futures[future] = chunk_i

            chunks_to_write = {}
            for future in as_completed(futures):
                chunks_to_write[futures[future]] = future.result()

        write_ordered_chunks(
            chunks_to_write,
            final_path=final_dir,
            append_dim="sample",
            engine="zarr",
            coord_dtypes=None,
            dtype=None,
        )

        final_da = xr.open_zarr(final_dir)
        subprocess.run(f"rm -rf {chunk_dir}", shell=True)
        update_dataset_config(output_dir=self.output_path, add_ds_region_dim={f"{region_dim}_prediction": region_dim},
                              change_region_dim=None)

        if region_dim == "query-region":
            self._region_prediction = final_da
        else:
            self._dmr_prediction = final_da
        return

    def predict(self, cpu, mask_cutoff=0.3, bw_bin_size=50):
        self._predict(region_dim="query-region", cpu=cpu, mask_cutoff=mask_cutoff)
        self._predict(region_dim="query-dmr", cpu=cpu, mask_cutoff=mask_cutoff)
        self.dump_bigwigs(cpu=cpu, mask_cutoff=mask_cutoff, bw_bin_size=bw_bin_size)
        return

    def _dump_sample(self, sample, mask_cutoff, bw_bin_size):
        # set dask scheduler to allow multiprocessing
        with dask.config.set(scheduler="sync"):
            dmr_pred = RegionDS.open(self.output_path, region_dim="query-dmr")
            region_pred = RegionDS.open(self.output_path, region_dim="query-region")

            # save DMR prediction proba
            dmr_bed_df = dmr_pred.get_bed(with_id=False)
            dmr_value = dmr_pred.get_feature(
                sample, dim="sample", da_name="query-dmr_prediction"
            )
            dmr_bed_df["score"] = dmr_value
            dmr_bed_df = dmr_bed_df[dmr_bed_df["score"] > mask_cutoff].copy()
            dmr_bed_df.sort_values(["chrom", "start"], inplace=True)
            dmr_bed_df.to_csv(
                f"{self.bigwig_dir}/{sample}_dmr_pred.bg",
                sep="\t",
                index=None,
                header=None,
            )

            # save region prediction proba
            region_bed_df = region_pred.get_bed(with_id=False)
            region_value = region_pred.get_feature(
                sample, dim="sample", da_name="query-region_prediction"
            )
            region_bed_df["score"] = region_value
            region_bed_df = region_bed_df[region_bed_df["score"] > mask_cutoff].copy()
            region_bed_df.sort_values(["chrom", "start"], inplace=True)
            region_bed_df.to_csv(
                f"{self.bigwig_dir}/{sample}_region_pred.bg",
                sep="\t",
                index=None,
                header=None,
            )

            bw_path = f"{self.bigwig_dir}/{sample}_reptile_score.bw"
            with pyBigWig.open(bw_path, "w") as bw:
                chrom_sizes = pd.read_csv(
                    self.chrom_size_path,
                    sep="\t",
                    index_col=0,
                    header=None,
                    squeeze=True,
                ).to_dict()
                bw.addHeader(
                    [(k, v) for k, v in pd.Series(chrom_sizes).sort_index().items()]
                )

                p = subprocess.run(
                    f"bedtools unionbedg -i "
                    f"{self.bigwig_dir}/{sample}_dmr_pred.bg "
                    f"{self.bigwig_dir}/{sample}_region_pred.bg",
                    shell=True,
                    check=True,
                    stdout=subprocess.PIPE,
                    encoding="utf8",
                )

                cur_bin = 0
                cur_scores = [0]
                for line in p.stdout.split("\n"):
                    if line == "":
                        continue
                    if line[-1] == ".":
                        # no score
                        continue

                    chrom, start, end, *scores = line.split("\t")
                    score = max(map(float, scores))
                    start_bin = int(start) // bw_bin_size
                    end_bin = int(end) // bw_bin_size + 1

                    for bin_id in range(start_bin, end_bin):
                        if bin_id > cur_bin:
                            # save previous bin
                            cur_pos = cur_bin * bw_bin_size
                            mean_score = sum(cur_scores) / len(cur_scores)
                            try:
                                bw.addEntries(
                                    chrom,
                                    [cur_pos],
                                    values=[mean_score],
                                    span=bw_bin_size,
                                )
                            except RuntimeError as e:
                                print(chrom, cur_pos, mean_score, bw_bin_size)
                                raise e

                            # init new bin
                            cur_bin = bin_id
                            cur_scores = [score]
                        elif bin_id == cur_bin:
                            # the same bin, take average
                            cur_scores.append(score)
                        else:
                            # no score, initial state
                            pass

                # final
                cur_pos = cur_bin * bw_bin_size
                mean_score = sum(cur_scores) / len(cur_scores)
                try:
                    bw.addEntries(
                        chrom, [cur_pos], values=[mean_score], span=bw_bin_size
                    )
                except RuntimeError as e:
                    print(chrom, cur_pos, mean_score, bw_bin_size)
                    raise e

            subprocess.run(
                f"rm -f {self.bigwig_dir}/{sample}_dmr_pred.bg "
                f"{self.bigwig_dir}/{sample}_region_pred.bg",
                shell=True,
            )
        return bw_path

    def dump_bigwigs(self, cpu, mask_cutoff, bw_bin_size):
        print(
            f"Save prediction results to bigwig files: \n"
            f"bin size {bw_bin_size}, minimum score {mask_cutoff}."
        )
        with ProcessPoolExecutor(cpu) as exe:
            futures = {}
            for sample in self.samples:
                f = exe.submit(
                    self._dump_sample,
                    sample=sample,
                    mask_cutoff=mask_cutoff,
                    bw_bin_size=bw_bin_size,
                )
                futures[f] = sample

            for f in as_completed(futures):
                sample = futures[f]
                print(f"{sample} result dump to: ", end="")
                bw_path = f.result()
                print(bw_path)
        return
