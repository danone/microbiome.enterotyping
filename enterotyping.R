#! /usr/bin/Rscript --vanilla

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  


# enterotyping R script based on DirichletMultinomial package
# usage example : Rscript --vanilla enterotyping.R -f otu_table_non_chimeric.biom -o enterotypes.txt


# set global option
library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="BIOM file name", metavar="character"),
	make_option(c("-o", "--out"), type="character", default="enterotypes.txt", 
              help="output file name [default= %default]", metavar="character"),
    make_option(c("-s", "--seed"), type="integer", default=444,
			  help="set seed number [default %default]", metavar="number")
          
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


# load libraries
library("biomformat")
library("magrittr")
library("DirichletMultinomial")
library("parallel")
library("reshape2")



# read arguments
input  = opt$file # for example from QIIME "otu_table_non_chimeric.biom"
output = opt$out # for example "enterotypes.txt"
seed   = opt$seed # for example 444, as default



# parse BIOM file and convert at genus levels
tax   = biomformat::read_biom(input) %>% biomformat::observation_metadata()
otu   = biomformat::read_biom(input) %>% biomformat::biom_data() %>% as.matrix()
genus = 
	otu %>% 
	apply(2,tapply,paste(
		tax[,1], 
		tax[,2], 
		tax[,3], 
		tax[,4], 
		tax[,5], 
		tax[,6]), 
	sum)



# enterotyping: fit a Dirichlet multinomial model
fit_genus_list = vector("list",5)


set.seed(seed); seeds=sample(1:1000, 5)

  for(i in 1:5) {

    set.seed(seeds[i])

    fit_genus <- mclapply(1:6, dmn, count=t(genus), verbose=FALSE, mc.cores=8)

    fit_genus_list[[i]] = fit_genus
    
    print(i)

  }


# collect Laplace score to find the best fit
lplc = vector("list",5)

for(i in 1:5) {

  lplc[[i]] <- sapply(fit_genus_list[[i]], function(x){attr(x,"goodnessOfFit")[["Laplace"]]})

}


# select the best number of cluster based on majority rule
best_genus_lplc = 
sapply(lplc, which.min) %>% table %>% which.max %>% names %>% as.integer


# assign enterotype id to each samples
enterotypes = 
fit_genus_list[[1]][[best_genus_lplc]] %>% 
mixture(assign=TRUE) %>% as.data.frame %>% set_colnames(c("Enterotypes_id"))




# write the output table
write.table(enterotypes, file=output, row.names=TRUE, sep="\t")








