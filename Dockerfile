FROM rocker/shiny

RUN apt-get update && apt-get install -y \
        libxml2-dev libssl-dev

RUN R -e "install.packages(c('data.table','httr','RCurl','shiny','shinyjs','shinyBS','DT','BiocManager', 'statmod'), repos='http://cran.cc.uoc.gr/mirrors/CRAN/')"
RUN R -e "install.packages(c('WebGestaltR','zip'))"
RUN R -e "BiocManager::install(c('GEOquery','limma','genefilter','annotate'))"
RUN R -e "BiocManager::install('oligo', configure.args='--disable-threading')"
RUN R -e "BiocManager::install('preprocessCore', configure.args='--disable-threading', force = TRUE)"

RUN apt-get install -y libharfbuzz-dev libfribidi-dev libudunits2-dev

RUN R -e "install.packages(c('bslib', 'htmltools', 'pillar', 'Rcpp', 'rlang', 'sass', 'textshaping', 'rematch2', 'colorspace', 'gridExtra', 'credentials', 'gitcreds', 'ini', 'farver', 'rstudioapi', 'brio', 'callr', 'desc', 'pkgload', 'praise', 'processx', 'ps', 'waldo', 'textshaping', 'labeling', 'munsell', 'RColorBrewer', 'viridisLite', 'downloader', 'influenceR', 'viridis', 'visNetwork', 'gtable', 'isoband', 'RcppCCTZ', 'zoo', 'RcppDate', 'diffobj', 'rex', 'gert', 'gh', 'rprojroot', 'thematic', 'testthat', 'markdown', 'Cairo', 'ragg', 'ape', 'igraph', 'igraphdata', 'rgl', 'scales', 'debugme', 'DiagrammeR', 'formattable', 'ggplot2', 'lubridate', 'nanotime', 'nycflights13', 'palmerpenguins', 'units', 'vdiffr', 'tinytest', 'inline', 'rbenchmark', 'pkgKitten', 'covr', 'usethis', 'RUnit'), repos='http://cran.cc.uoc.gr/mirrors/CRAN/')"

RUN R -e "BiocManager::install('graph')"

COPY /DExplore /home/rstudio/
COPY /DExplore /srv/shiny-server/dexplore/

RUN chown -R shiny /srv/shiny-server/dexplore/
RUN chown -R shiny /usr/local/lib/R/site-library

RUN chmod -R +w /srv/shiny-server/dexplore/
RUN chmod -R +w /usr/local/lib/R/site-library

EXPOSE 3838
EXPOSE 8787