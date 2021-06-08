FROM conda/miniconda3-centos7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/chip-seq.git"
LABEL gitCommit="201185c7f380bcdc31f6d4b3cd44f7086fcd238a / devel"

# real path from baseDir: /home/nservant/Apps/geniac/build/workDir/recipes/conda/bowtie2.yml
ADD bowtie2.yml /opt/bowtie2.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/bowtie2.yml \
&& echo "source activate bowtie2" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/bowtie2/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
