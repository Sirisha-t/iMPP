FROM ubuntu AS base

MAINTAINER sirisha

WORKDIR /project/impp

##build layer
FROM base AS build

##update and install dependencies
RUN apt-get update -qy
RUN apt-get upgrade -qy
RUN apt-get install build-essential -qy
RUN apt-get install perl
RUN apt-get install python -qy
#RUN apt-get install python3 -qy
RUN apt-get install libboost-dev -qy
RUN apt-get install zlib1g-dev -qy


##runtime layer
FROM base AS run

COPY --from=build /usr/local /usr/local

##copy directories
COPY --chown=docker:docker  install.sh /project/impp/install.sh
COPY --chown=docker:docker README /project/impp/README
COPY --chown=docker:docker /lib /project/impp/lib
COPY --chown=docker:docker /src /project/impp/src
COPY --chown=docker:docker /utils /project/impp/utils
COPY --chown=docker:docker /params /project/impp/params
COPY --chown=docker:docker impp_run.pl /project/impp/impp_run.pl


RUN echo 'we are running some # of cool things'
##create dir for executables
RUN mkdir -p /project/impp/bin
##run install script
RUN bash /project/impp/install.sh


