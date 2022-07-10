def install_r_package(name):
    from rpy2.robjects.vectors import StrVector
    from rpy2.robjects.packages import importr, isinstalled

    if not isinstalled(name):
        utils = importr("utils")
        utils.chooseCRANmirror(ind=1)
        utils.install_packages(StrVector([name]))


def install_github_r_package(github_name):
    from rpy2.robjects.packages import importr, isinstalled

    install_r_package("devtools")
    devtools = importr("devtools")

    if not isinstalled(github_name.split('/')[-1]):
        devtools.install_github(github_name)
    return

