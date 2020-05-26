.. currentmodule:: MAPT


Managing Sequence Data from Silva
----------------------------------------

.. autoclass:: MAPT.silva_manager
    :members:

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Mapping Kmers 
----------------------------------------

.. autoclass:: MAPT.k_mapper
    :members:

.. toctree::
   :maxdepth: 2
   :caption: Contents:

R code for k-mer mapping figures in docs
------------------------------------------
Feel free to use this R code to generate plots showing k-mer mappings!
~~~
PNA <- read.csv('path/to/pna/results/after/mapping/pna')

library('ggplot2')
axis_text_size = 14

org_name <- "your organism name"

plt<-ggplot(PNA,aes(x=index,y=absolute.match))+
    geom_line()+
    geom_line(aes(y=unique.match),color='blue')+
    geom_line(aes(y=PNA.mapping),color = 'red')+
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(.x)))+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
       axis.line = element_line(color="black"),
       axis.text = element_text(size=axis_text_size),
       axis.title=element_text(size=axis_text_size))+
    labs(y="Log K-mer Matches",x=org_name)
~~~
