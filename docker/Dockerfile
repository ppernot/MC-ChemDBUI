FROM ppernot1/shiny_base

RUN apt-get install -y --no-install-recommends libharfbuzz-dev libfribidi-dev

# renv
ENV RENV_VERSION 0.15.5
RUN Rscript -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN Rscript -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

RUN mkdir MC-ChemDBUI
WORKDIR /MC-ChemDBUI
COPY .Rprofile .
#ADD renv.tar.gz . # Pb unresolved shared libs...
COPY renv/activate.R ./renv/activate.R
COPY renv.lock ./renv.lock
RUN R -e 'renv::consent(TRUE)' \
 && R -e 'renv::restore()'

## database files
WORKDIR /
RUN mkdir ChemDBPublic
RUN mkdir MC-ChemDB
#ADD ChemDB.tar.gz /  ### Now linked to user's disk

## Sampler files
WORKDIR /MC-ChemDBUI
COPY Sampler ./Sampler/
WORKDIR /MC-ChemDBUI/Sampler/Ranlib/src
RUN ./comp
WORKDIR /MC-ChemDBUI/Sampler
RUN ./compile

## app files
WORKDIR /MC-ChemDBUI
COPY global.R .
COPY ui.R .
COPY server.R .
COPY ui_files ./ui_files/
COPY server_files ./server_files/
COPY R ./R/
COPY data ./data/

## change user from root ### BREAKS renv links !!!
#RUN addgroup --system app \
#    && adduser --system --ingroup app app
#RUN chown app:app -R /UncVal
#USER app

## expose port and run app
EXPOSE 3820
CMD ["R", "-e", "shiny::runApp(host = '0.0.0.0', port = 3820)"]
