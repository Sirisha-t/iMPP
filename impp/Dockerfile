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
RUN apt-get install python3 -qy
RUN apt-get install libboost-dev -qy
RUN apt-get install zlib1g-dev -qy


##runtime layer
FROM base AS run

COPY --from=build /usr/local /usr/local

##copy directories
COPY install.sh /project/impp/install.sh
COPY README /project/impp/README
COPY /lib /project/impp/lib
COPY /src /project/impp/src
COPY /utils /project/impp/utils
COPY /params /project/impp/params
COPY  impp_run.pl /project/impp/impp_run.pl


##create dir for executables
RUN mkdir -p /project/impp/bin
##run install script
RUN bash /project/impp/install.sh


