FROM bioconductor/bioconductor_docker:RELEASE_3_23

ARG repo_name=daisymut/sammyR
ARG tag_name

RUN R -e 'install.packages(c("data.table","patchwork","dendextend","effsize","digest","remotes"))'
RUN R -e 'BiocManager::install(c("GenomicRanges","rtracklayer","S4Vectors","GenomeInfoDb","Gviz","GenomicFeatures","limma","IRanges","BiocGenerics","txdbmaker"), update=FALSE, ask=FALSE)'
RUN R -e 'remotes::install_github("CSOgroup/CALDER2", upgrade="never")'
RUN R -e "remotes::install_github('${repo_name}@${tag_name}', upgrade='never')"

LABEL org.opencontainers.image.source=https://github.com/daisymut/sammyR
