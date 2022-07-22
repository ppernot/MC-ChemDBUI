[![DOI](https://zenodo.org/badge/515063130.svg)](https://zenodo.org/badge/latestdoi/515063130)

# MC-ChemDBUI

UI to manage [MC-ChemDB](https://github.com/ppernot/MC-ChemDB).

## News

* `2022/07/21` First docker version to handle neutral reactions with three main functions:

    1. Open and edit DB files; create new DB version

    2. Parse the DB files to check the consistency (names, masses, mass balance...)
    
    3. Generate MC samples and plot resulting reaction rates as a function of temperature
    and density.
    
## To Be Done...

* Add possibility to Import/Export MC-ChemDB (to and from github ?)

* Manage ion-molecule reactions

* Manage photo-processes


## Docker container

The [mc-chemdbui](https://hub.docker.com/repository/docker/ppernot1/mc-chemdbui)
Docker container has all elements pre-installed, including [MC-ChemDB](https://github.com/ppernot/MC-ChemDB).

To run the container:

0. Install [Docker](https://www.docker.com/products/docker-desktop)

1. Type the following command in a terminal
```
docker run -d -p 3820:3820 --mount type=bind,source="$(pwd)"/../ChemDBPublic,target=/ChemDBPublic --name mc-chemdbui ppernot1/mc-chemdbui
```
where `ChemDBPublic` is a local directory containing the processed MC samples

2. Access MC-ChemDBUI at `http://localhost:3820` in your favorite browser

3. When finished
```
docker kill mc-chemdbui
```
__Warning__: all modifs to the database will be lost. To keep them,
you need to create a new image of the modified container...
```
docker commit mc-chemdbui mc-chemdbui:vX.X
```
where `X.X` is the new tag. 
You should then refer to this specific version... 

4. For further sessions
```
docker restart mc-chemdbui
```

4. To cleanup
```
docker remove -v mc-chemdbui
```
