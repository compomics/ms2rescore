FROM debian:stable-slim as percolator
RUN apt-get update && apt-get install -y --no-install-recommends wget tar gzip ca-certificates
WORKDIR /tmp
ARG PERCOLATOR_VERSION=3-05
RUN wget https://github.com/percolator/percolator/releases/download/rel-${PERCOLATOR_VERSION}/ubuntu.tar.gz && \
    tar -xzf ubuntu.tar.gz

FROM python:3.8-buster
ARG PERCOLATOR_VERSION=3-05
COPY --from=percolator /tmp/percolator-noxml-v${PERCOLATOR_VERSION}-linux-amd64.deb /tmp
COPY --from=percolator /tmp/percolator-converters-v${PERCOLATOR_VERSION}-linux-amd64.deb /tmp
RUN dpkg -i /tmp/percolator-noxml-v${PERCOLATOR_VERSION}-linux-amd64.deb && \
    dpkg -i /tmp/percolator-converters-v${PERCOLATOR_VERSION}-linux-amd64.deb && \
    rm /tmp/percolator-noxml-v${PERCOLATOR_VERSION}-linux-amd64.deb && \
    rm /tmp/percolator-converters-v${PERCOLATOR_VERSION}-linux-amd64.deb
COPY . .
RUN pip install -r requirements.txt
EXPOSE 8501
CMD ["streamlit", "run", "ms2rescore_streamlit.py"]
