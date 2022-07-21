[![DOI](https://zenodo.org/badge/515063130.svg)](https://zenodo.org/badge/latestdoi/515063130)

# MC-ChemDBUI

UI for MC-ChemDB

## Docker container

The [mc-chemdbui](https://hub.docker.com/repository/docker/ppernot1/mc-chemdbui)
Docker container has all elements preinstalled.

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
you need to create a new image of the modified container... (TBD)

4. For further sessions
```
docker restart mc-chemdbui
```

4. To cleanup
```
docker remove -v mc-chemdbui
```
