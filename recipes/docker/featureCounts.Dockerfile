FROM conda/miniconda3-centos7

LABEL gitUrl="ssh://git@gitlab.curie.fr:2222/data-analysis/chip-seq.git"
LABEL gitCommit="201185c7f380bcdc31f6d4b3cd44f7086fcd238a / devel"

# real path from baseDir: /home/nservant/Apps/geniac/build/workDir/recipes/conda/featureCounts.yml
ADD featureCounts.yml /opt/featureCounts.yml

RUN yum install -y which   \
&& yum clean all \
&& conda env create -f /opt/featureCounts.yml \
&& echo "source activate featureCounts" > ~/.bashrc \
&& conda clean -a


ENV PATH /usr/local/envs/featureCounts/bin:$PATH

ENV LC_ALL en_US.utf-8
ENV LANG en_US.utf-8
