FROM ppernot1/shiny_base

RUN apt-get install -y --no-install-recommends libharfbuzz-dev libfribidi-dev

## renv
ENV RENV_VERSION 0.15.5
RUN Rscript -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN Rscript -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

RUN mkdir MC-ChemDBUI
COPY .Rprofile /MC-ChemDBUI/
COPY renv/activate.R /MC-ChemDBUI/renv/activate.R
COPY renv.lock /MC-ChemDBUI/renv.lock
WORKDIR /MC-ChemDBUI
RUN R -e 'renv::consent(TRUE)' \
 && R -e 'renv::restore()'

## database files
WORKDIR /
RUN mkdir ChemDBPublic
ADD ChemDB.tar.gz /

## app files
COPY global.R /MC-ChemDBUI/
COPY ui.R /MC-ChemDBUI/
COPY server.R /MC-ChemDBUI/
COPY ui_files /MC-ChemDBUI/ui_files/
COPY server_files /MC-ChemDBUI/server_files/
COPY R /MC-ChemDBUI/R/
COPY data /MC-ChemDBUI/data/
WORKDIR /MC-ChemDBUI


## change user from root ### BREAKS renv links !!!
#RUN addgroup --system app \
#    && adduser --system --ingroup app app
#RUN chown app:app -R /UncVal
#USER app

## expose port and run app
EXPOSE 3820
CMD ["R", "-e", "shiny::runApp(host = '0.0.0.0', port = 3820)"]
