FROM centos:7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/chip-seq.git"
LABEL gitCommit="201185c7f380bcdc31f6d4b3cd44f7086fcd238a / devel"

RUN yum install -y which \
&& yum clean all

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
