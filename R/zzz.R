.onAttach <- function(libname, pkgname)
{
    packageStartupMessage( paste0(
        "Because of space limitations, TFEA.ChIPs internal database only includes ",
        "ChIP-seq experiments from cell types in ENCODE's tiers 1, 2, and 2.5. \n",
        "To download the full ReMap2022 database, as well as other ready-to-use databases, ",
        "visit https://github.com/LauraPS1/TFEA.ChIP_downloads"
    )
    )
}
