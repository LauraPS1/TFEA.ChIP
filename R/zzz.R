.onAttach <- function(libname, pkgname)
{
    packageStartupMessage( paste0(
        "Because of space limitations, TFEA.ChIPs internal database only includes ",
        "ChIP-seq experiments from the ENCODE project. \n",
        "To download the full ReMap database, as well as other ready-to-use databases, ",
        "visit https://github.com/LauraPS1/TFEA.ChIP_downloads"
    )
    )
}
