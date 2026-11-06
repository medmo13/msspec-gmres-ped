# Get the base Python image
FROM alpine:edge AS builder

# Variables
ARG branch="devel"
ARG login="" password=""
ARG folder=/opt/msspec user=msspec

# Install system dependencies
# tools
RUN apk add bash git make gfortran python3 py3-numpy-f2py 
# headers
RUN apk add python3-dev lapack-dev musl-dev hdf5-dev cairo-dev
# python packages
RUN apk add py3-virtualenv py3-pip py3-numpy-dev py3-h5py py3-lxml py3-matplotlib \
            py3-numpy py3-pandas py3-cairo py3-scipy py3-setuptools_scm \
            py3-terminaltables ipython
RUN apk add --no-cache -X http://dl-cdn.alpinelinux.org/alpine/edge/community py3-wxpython
#RUN pip install ase pint terminaltables ipython
# for GUI
RUN apk add ttf-droid adwaita-icon-theme
RUN apk add build-base

# Fetch the code
RUN mkdir -p ${folder}/code
WORKDIR ${folder}/code
RUN git clone --branch ${branch} https://${login}:${password}@git.ipr.univ-rennes1.fr/epsi/msspec_python3.git .

RUN virtualenv --system-site-packages ${folder}/.local/src/msspec_venv
RUN make        pybinding PYTHON=python3 VENV_PATH=${folder}/.local/src/msspec_venv VERBOSE=1
RUN make -C src sdist     PYTHON=python3 VENV_PATH=${folder}/.local/src/msspec_venv VERBOSE=1
RUN make -C src frontend  PYTHON=python3 VENV_PATH=${folder}/.local/src/msspec_venv VERBOSE=1
RUN source ${folder}/.local/src/msspec_venv/bin/activate && pip install src/dist/msspec*tar.gz



# Build
#RUN make pybinding NO_VENV=1 PYTHON=python3 VERBOSE=1
#RUN make -C src sdist PYTHON=python3 NO_VENV=1 VENV_PATH=${folder}/.local/src/msspec_venv 
#&& \
#    pip install src/dist/msspec*tar.gz

# Add a non-privileged user
#RUN adduser -D -s /bin/bash -h ${folder} ${user}

# Set the working directory in the container
#USER ${user}

#RUN virtualenv --system-site-packages ${folder}/.local/src/msspec_venv
#RUN source ${folder}/.local/src/msspec_venv/bin/activate && pip install src/dist/msspec*.tar.gz
#RUN make -C src frontend PYTHON=python3 NO_VENV=1 VENV_PATH=${folder}/.local/src/msspec_venv



FROM alpine:edge
# Variables
ARG folder=/opt/msspec user=msspec
# Install system dependencies
RUN apk add --no-cache -X http://dl-cdn.alpinelinux.org/alpine/edge/community \
#  hdf5-hl cairo openblas lapack libxml2 libxslt libzlf wxwidgets-gtk3 openjpeg libimagequant \
   nano \
   py3-virtualenv \
   lapack \
   bash \
#  git \
#  make \
#  gfortran \
   python3 \
#   ttf-droid \
   ttf-liberation \
   adwaita-xfce-icon-theme \
#  python3-dev \
#  lapack-dev \
#  musl-dev \
#  py3-virtualenv \
   py3-pip \
#  py3-numpy-dev \
   py3-h5py \
   py3-lxml \
   py3-matplotlib \
   py3-numpy \
   py3-pandas \
   py3-cairo \
   py3-scipy \
   py3-setuptools_scm \
   py3-wxpython \
   py3-terminaltables \
   py3-bayesian-optimization \
   # Add a non-privileged user
   && adduser -D -s /bin/bash -h ${folder} ${user}

# Set the working directory in the container
USER ${user}
WORKDIR ${folder}
# Install msspec
#COPY --from=builder ${folder}/.local ${folder}/.local
#COPY --from=builder /usr/lib/python3.10/site-packages /usr/lib/python3.10/site-packages

COPY --from=builder ${folder}/code/src/dist/msspec*tar.gz msspec.tar.gz
RUN virtualenv --system-site-packages .local/src/msspec_venv && \
    . .local/src/msspec_venv/bin/activate && \
    pip install msspec.tar.gz && \
    pip install ipython && \
    pip cache purge && \
    rm -f msspec.tar.gz && \
    mkdir -p .local/bin

COPY --from=builder /root/.local/bin/msspec .local/bin/msspec
ENV PATH=${folder}/.local/bin:$PATH

# Run the msspec frontend command on startup
ENTRYPOINT ["msspec"]
