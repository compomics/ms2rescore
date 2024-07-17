FROM python:3.11

# ARG DEBIAN_FRONTEND=noninteractive

LABEL name="ms2rescore"

# ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ms2rescore

ADD pyproject.toml /ms2rescore/pyproject.toml
ADD LICENSE /ms2rescore/LICENSE
ADD README.md /ms2rescore/README.md
ADD MANIFEST.in /ms2rescore/MANIFEST.in
ADD ms2rescore /ms2rescore/ms2rescore

RUN apt-get update \
    && apt install -y procps \
    && pip install /ms2rescore --only-binary :all:

ENTRYPOINT [""]
