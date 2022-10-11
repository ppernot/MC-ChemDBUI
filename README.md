[![DOI](https://zenodo.org/badge/515063130.svg)](https://zenodo.org/badge/latestdoi/515063130)

# MC-ChemDBUI

UI to manage [MC-ChemDB](https://github.com/ppernot/MC-ChemDB).

## News

* `2022/07/21` First docker version (v0.1) to handle neutral reactions with three main functions:

    1. Open and edit DB files; create new DB version

    2. Parse the DB files to check the consistency (names, masses, mass balance...)
    
    3. Generate MC samples and plot resulting reaction rates as a function of temperature
    and density.

* `2022/07/29` v0.2

    1. Add ion/molecule reactions (edit/parse, sample). No `Save` feature yet...

* `2022/09/13` v0.3

    1. Neutrals and Ions are now functional, with similar interfaces enabling to edit and save `MC-ChemDB` database and to generate samples in `ChemDBPublic`.
    
    2. Uses a new version of [MC-ChemDB](https://github.com/ppernot/MC-ChemDB) with a single .csv file for neutrals and another one for ions. This enables an easier management of [MC-ChemDB](https://github.com/ppernot/MC-ChemDB) with GitHub.
    
* `2022/10/11` v0.4

    1. Implements the management of photo-processes.
    
    
## Documentation

See the [User's Manual](https://raw.githubusercontent.com/ppernot/MC-ChemDBUI/main/docs/manual.pdf)

## Docker container

The [mc-chemdbui](https://hub.docker.com/repository/docker/ppernot1/mc-chemdbui)
Docker container has all elements pre-installed.

To run the container:

0. Install [Docker](https://www.docker.com/products/docker-desktop)

1. Load the latest release of `MC-ChemDB` from the GitHub repository https://github.com/ppernot/MC-ChemDB and unpack it 

2. If it does not exist, create a `ChemDBPublic` directory 

3. Download the latest docker container for `mc-chemdbui` from DockerHub 
```
docker pull ppernot1/mc-chemdbui
``` 

4. Run the docker container with source links for MC-ChemDB and ChemDBPublic pointing to *your* directories 
```
docker run -d -p 3820:3820\
	--mount type=bind,source=path_to_my_ChemDBPublic,target=/ChemDBPublic\
	--mount type=bind,source=path_to_my_MC-ChemDB,target=/MC-ChemDB\
	--name mc-chemdbui ppernot1/mc-chemdbui 
```

5. Access http://localhost:3820 in your favorite browser

6. When finished
```
docker kill mc-chemdbui
```

7. For further sessions
```
docker restart mc-chemdbui
```

8. To cleanup
```
docker remove -v mc-chemdbui
```
