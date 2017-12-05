#Variant analysis on a single sample 
#Pavitra Roychoudhury, Aug 2017

rm(list=ls()); 
sessionInfo();
library(ShortRead); library(parallel); 
ncores<-detectCores();
# ncores<-8;

#Get args from command line 
args<-(commandArgs(TRUE));
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
    print(args[[i]])
  }
}
#For testing
#paired<-TRUE; #from command line
#tgt_region<-'ATAAACTCACACACGGCGTCCTGG'; #should come from command line
##s1 and s2 for paired end reads, or just s1 for single-end
#s1<-'Martine-0417-1_S351_L001_R1_001.fastq.gz'
#s2<-'Martine1_S327_L001_R2_001.fastq.gz'
#refseq_fname<-'./hsv1m5.fasta'

# paired<-FALSE; #from command line
# tgt_region<-'GGCCCTTTGACGCCGAGACCAGACGGGT'; #should come from command line
# s1<-'martine-5_S221_L001_R1_001.fastq.gz'
# refseq_fname<-'./hsv_crispr_amp.fasta'


#Files, directories, target site
bamfdir<-'./mapped_reads/'; 
fastq_dir<-'./fastq_files/'; 
output_dir<-'./results/';
sampname<-strsplit(s1,'_R1_001.fastq.gz')[[1]][1];

# Run QA on fastq files and extract read counts 
# report(qaSummary,dest=paste(output_dir,s1,'_QAReport',sep=''),type='html')
if(!paired){
  qaSummary<-qa(fastq_dir,pattern=s1,type='fastq');
  n_rawreads<-qaSummary[['readCounts']]$read;
}else{
  qaSummary1<-qa(fastq_dir,pattern=s1,type='fastq');
  qaSummary2<-qa(fastq_dir,pattern=s2,type='fastq');
  n_rawreads<-qaSummary1[['readCounts']]$read+qaSummary2[['readCounts']]$read;
}


## ----import_data---------------------------------------------------------
ref_seq<-readFasta(refseq_fname);
bamfname<-grep(sampname,list.files(bamfdir,'*.sorted.bam$'),value=T);
baifname<-indexBam(paste(bamfdir,bamfname,sep='')); #Make an index file
params<-ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE),
                     what=c('qname','rname','strand','pos','qwidth','mapq','cigar','seq','qual'));
reads<-scanBam(file=paste(bamfdir,bamfname,sep=''),
               index=paste(bamfdir,bamfname,'.bai',sep=''),
               param=params); #import bam file

#locate target sites in ref seq and get start and end positions
tgt_pos<-vmatchPattern(tgt_region,sread(ref_seq))[[1]];  

## ----read_preview--------------------------------------------------------
#     names(reads[[1]]); #this lists all the elements in the bam file
#     head(reads[[1]]$qname); #read IDs
#     head(reads[[1]]$rname); #name of reference seq
#     head(reads[[1]]$strand); #+ or - strand
#     plot(reads[[1]]$pos,main='POS field in BAM file',ylab='pos'); #leftmost position in reference of first matching base
#     hist(reads[[1]]$pos,xlab='First matching position')
#     head(reads[[1]]$qwidth); #read lengths
#     head(reads[[1]]$mapq); #map quality 
#     head(reads[[1]]$cigar); #the alignment
#     reads[[1]]$seq; #sequences in the bam file
#     reads[[1]]$qual; #quality scores
#     hist(alphabetScore(PhredQuality(reads[[1]]$qual))/width(PhredQuality(reads[[1]]$qual)),
#          main='Avg read quality',xlab='')

## ----qual_filter---------------------------------------------------------
goodavgqual<-srFilter(function(x){
  alphabetScore(x$quals)/x$widths>30;
},name='GoodAvgQual')
good_qual<-goodavgqual(list(quals=reads[[1]]$qual,widths=reads[[1]]$qwidth));
good_qual;

## ----cover_target--------------------------------------------------------
covertgt<-srFilter(function(x){
  x$start<=start(tgt_pos)&
    x$end>=end(tgt_pos)
},name='CoverTargetRegion')
cover_tgt<-covertgt(list(start=reads[[1]]$pos,end=(reads[[1]]$pos+reads[[1]]$qwidth-1)));
cover_tgt;

## ----all_filters---------------------------------------------------------
filter_inds<-as.logical(good_qual)&as.logical(cover_tgt);
sum(filter_inds);

## ----identify reads with variants in tgt----------------------------------------
varread_inds<-filter_inds&rowSums(cigarOpTable(reads[[1]]$cigar)[,c('I','D')])>0;
sum(varread_inds);

## ----freq_table----------------------------------------------------------
temp_table<-tables(reads[[1]]$seq[varread_inds],
                   n=length(unique(reads[[1]]$seq[varread_inds])));
freq_table<-data.frame(seq=names(temp_table$top),
                       length=nchar(names(temp_table$top)),
                       freq=unname(temp_table$top),stringsAsFactors=F);
rm(temp_table);
head(freq_table);
tail(freq_table);

## ----also build and save a pileup----------------------------------------
pileupdf<-pileup(file=paste(bamfdir,bamfname,sep=''),
                 index=paste(bamfdir,bamfname,'.bai',sep=''),
                 pileupParam=PileupParam(max_depth=sum(filter_inds),min_base_quality=30,min_mapq=0,
                                         min_nucleotide_depth=1,min_minor_allele_depth=0,
                                         distinguish_strands=T,distinguish_nucleotides=T,
                                         ignore_query_Ns=T,include_deletions=T,include_insertions=T));

if(!dir.exists('./pileups/')) dir.create('./pileups/');
write.csv(pileupdf,file=paste('./pileups/',sampname,'_pileup.csv',sep='')); rm(pileupdf)

## ----variant_analysis----------------------------------------------------
var_table<-freq_table;
ops<-setdiff(CIGAR_OPS, c('P','M')); #ignore padding and matches
inds<-unlist(mclapply(c(1:nrow(var_table)),function(x)
  match(var_table$seq[x],reads[[1]]$seq[varread_inds]),mc.cores=ncores)); #takes a while
cigs_ref<-cigarRangesAlongReferenceSpace(reads[[1]]$cigar[varread_inds][inds],ops=ops,with.ops=T,
                                         reduce.ranges=T,drop.empty.ranges=F,
                                         pos=reads[[1]]$pos[varread_inds][inds]);
cigs_query<-cigarRangesAlongQuerySpace(reads[[1]]$cigar[varread_inds][inds],ops=ops,with.ops=T,
                                       reduce.ranges=T,drop.empty.ranges=F);
cigs_pair<-cigarRangesAlongPairwiseSpace(reads[[1]]$cigar[varread_inds][inds],ops=ops,with.ops=T,
                                         reduce.ranges=T,drop.empty.ranges=F);
cigs_query_offset<-mclapply(c(1:length(cigs_pair)),
                            function(i)return(IRanges(start=start(cigs_pair)[[i]]+
                                                        reads[[1]]$pos[varread_inds][inds[i]]-1,
                                                      end=end(cigs_pair[[i]])+
                                                        reads[[1]]$pos[varread_inds][inds[i]]-1,
                                                      width=width(cigs_pair[[i]]),
                                                      names=names(cigs_pair[[i]]))),mc.cores=ncores);
labels_cigs<-mclapply(cigs_ref,function(x)names(x),mc.cores=ncores);
tmp_labels_indels<-mclapply(labels_cigs,function(x)unlist(mclapply(x,function(y)
  strsplit(y,'')[[1]][1],mc.cores=ncores)),mc.cores=ncores); #reduces longer indel labels to 1 char
rm(labels_cigs);
var_table$start_pos<-reads[[1]]$pos[varread_inds][inds];

#Count indels
var_table$num_insertions<-unlist(mclapply(tmp_labels_indels,function(x)sum(x=='I'),mc.cores=ncores));
var_table$num_deletions<-unlist(mclapply(tmp_labels_indels,function(x)sum(x=='D'),mc.cores=ncores));

#Length and positions of deletions
var_table$del_frameshift<-var_table$del_in_target<-
  var_table$len_deletions<-var_table$pos_deletions<-'';
pos_del<-mclapply(tmp_labels_indels,function(x)which(x=='D'),mc.cores=ncores); 
locs<-which(unlist(mclapply(pos_del,function(x)length(x)>0,mc.cores=ncores))); 
if(length(locs)>0){
  tmp<-mclapply(locs,function(x)start(cigs_ref)[[x]][pos_del[[x]]],mc.cores=ncores); #takes a while
  var_table$pos_deletions[locs]<-unlist(mclapply(tmp,function(x)paste(x,collapse=';'),mc.cores=ncores));
  tmp<-mclapply(locs,function(x)width(cigs_ref)[[x]][pos_del[[x]]],mc.cores=ncores); #takes a while
  var_table$len_deletions[locs]<-unlist(mclapply(tmp,function(x)paste(x,collapse=';'),mc.cores=ncores));
  tmp_del<-mclapply(locs,function(x)cigs_ref[[x]][pos_del[[x]]],mc.cores=ncores); 
  hits_del<-mclapply(tmp_del,function(del)countOverlaps(del,tgt_pos),mc.cores=ncores);
  var_table$del_in_target[locs]<-unlist(mclapply(hits_del,function(x)
    paste(as.logical(x),collapse=';'),mc.cores=ncores));
  var_table$del_frameshift[locs]<-unlist(mclapply(var_table$len_deletions[locs],function(x)
    paste(as.numeric(strsplit(x,';')[[1]])%%3!=0,collapse=';'),mc.cores=ncores));
  rm(tmp,locs,tmp_del,hits_del); gc();
}

#Length and positions of insertions
var_table$ins_frameshift<-var_table$ins_in_target<-
  var_table$len_insertions<-var_table$pos_insertions<-'';
pos_ins<-mclapply(tmp_labels_indels,function(x)which(x=='I'),mc.cores=ncores);
locs<-which(unlist(mclapply(pos_ins,function(x)length(x)>0,mc.cores=ncores))); 
if(length(locs)>0){
  tmp<-mclapply(locs,function(x)start(cigs_ref)[[x]][pos_ins[[x]]],mc.cores=ncores); #takes a while
  var_table$pos_insertions[locs]<-unlist(mclapply(tmp,function(x)paste(x,collapse=';'),mc.cores=ncores)); 
  tmp<-mclapply(locs,function(x)width(cigs_query)[[x]][pos_ins[[x]]],mc.cores=ncores);
  var_table$len_insertions[locs]<-unlist(mclapply(tmp,function(x)paste(x,collapse=';'),mc.cores=ncores));
  tmp_ins<-mclapply(locs,function(x)cigs_query_offset[[x]][pos_ins[[x]]],mc.cores=ncores); #takes a while
  hits_ins<-mclapply(tmp_ins,function(ins)countOverlaps(ins,tgt_pos),mc.cores=ncores);
  var_table$ins_in_target[locs]<-unlist(mclapply(hits_ins,function(x)
    paste(as.logical(x),collapse=';'),mc.cores=ncores));
  var_table$ins_frameshift[locs]<-unlist(mclapply(var_table$len_insertions[locs],function(x)
    paste(as.numeric(strsplit(x,';')[[1]])%%3!=0,collapse=';'),mc.cores=ncores));
  rm(tmp,tmp_labels_indels);gc();
}

#Sequence and quality of insertions
var_table$errprob_insertions<-var_table$phredqual_insertions<-
  var_table$char_insertions<-'';
if(length(locs)>0){
  tmp_ranges<-mclapply(locs,function(x)cigs_query[[x]][pos_ins[[x]]],mc.cores=ncores);
  tmp_seqs<-reads[[1]]$seq[varread_inds][inds][locs];
  var_table$char_insertions[locs]<-unlist(mclapply(c(1:length(locs)),function(i)
    paste(mclapply(tmp_ranges[[i]],function(x)as.character(tmp_seqs[[i]][x]),mc.cores=ncores),
          collapse=';'),mc.cores=ncores));
  rm(tmp_seqs);
  tmp_quals<-reads[[1]]$qual[varread_inds][inds][locs];
  tmp_qual<-unlist(mclapply(c(1:length(locs)),function(i)mclapply(tmp_ranges[[i]],function(x)
    PhredQuality(tmp_quals[[i]][x]),mc.cores=ncores),mc.cores=ncores)); rm(tmp_quals);
  var_table$phredqual_insertions[locs]<-unlist(mclapply(c(1:length(locs)),function(i)
    paste(unlist(mclapply(tmp_qual[[i]],function(x)
      round(alphabetScore(PhredQuality(x))/width(PhredQuality(x)),1),mc.cores=ncores)),
      collapse=';'),mc.cores=ncores));
  rm(tmp_qual); tmp_phred<-var_table$phredqual_insertions[locs];
  var_table$errprob_insertions[locs]<-unlist(mclapply(tmp_phred,function(x)
    paste(round(10^(-as.numeric(strsplit(x,';')[[1]])/10),digits=4),collapse=';'),mc.cores=ncores));
  rm(tmp_ranges,tmp_phred,locs); gc()
}

var_table$name<-reads[[1]]$qname[varread_inds][inds];
head(var_table);
write.csv(var_table,file=paste(output_dir,sampname,'_variantanalysis.csv',sep=''),row.names=F)

## ----summarize-----------------------------------------------------------
var_table<-read.csv(file=paste(output_dir,sampname,'_variantanalysis.csv',sep=''),
                    stringsAsFactors=F);

## need to check the insertion numbers, they're not adding up
var_summary<-data.frame(
  file=sampname,
  pairedend=paired,
  tgt=tgt_region,
  orig_reads=n_rawreads,
  mapped_reads=length(reads[[1]]$qname),
  good_qual=stats(good_qual)$Passing,
  cover_tgt=stats(cover_tgt)$Passing,
  tot_after_filters=sum(filter_inds),
  reads_with_vars=sum(varread_inds),
  totreads_ins=sum(var_table$freq[var_table$num_insertions>0]),
  totreads_del=sum(var_table$freq[var_table$num_deletions>0]),
  totreads_var=sum(var_table$freq[var_table$num_insertions>0|var_table$num_deletions>0]),
  totreads_del_tgt=sum(var_table$freq[var_table$num_deletions>0&
                                        grepl('TRUE',var_table$del_in_target)]),
  totreads_ins_tgt=sum(var_table$freq[var_table$num_insertions>0&
                                        grepl('TRUE',var_table$ins_in_target)]),
  totreads_var_tgt=sum(var_table$freq[(var_table$num_insertions>0|var_table$num_deletions>0)&
                                        (grepl('TRUE',var_table$ins_in_target)|
                                           grepl('TRUE',var_table$del_in_target))]),
  totreads_ins_tgt_inframe=sum(var_table$freq[grepl('TRUE',var_table$ins_in_target)&
                                                grepl('FALSE',var_table$ins_frameshift)]),
  totreads_ins_tgt_frameshift=sum(var_table$freq[grepl('TRUE',var_table$ins_in_target)&
                                                   grepl('TRUE',var_table$ins_frameshift)]),
  totreads_del_tgt_inframe=sum(var_table$freq[grepl('TRUE',var_table$del_in_target)&
                                                grepl('FALSE',var_table$del_frameshift)]),
  totreads_del_tgt_frameshift=sum(var_table$freq[grepl('TRUE',var_table$del_in_target)&
                                                   grepl('TRUE',var_table$del_frameshift)])
	totreads_novariants_in_tgt=sum(var_table$freq[grepl('FALSE',var_table$del_in_target)&
																									grepl('FALSE',var_table$del_frameshift)]));
var_summary;
write.csv(var_summary,file=paste(output_dir,sampname,'_variantsummary.csv',sep=''),row.names=F);

## ----export_reads--------------------------------------------------------

#Write to fasta
# reads_out<-DNAStringSet(var_table$seq);
# names(reads_out)<-var_table$name;
# writeXStringSet(reads_out,filepath=paste(output_dir,sampname,'_analyzedreads.fa',sep=''));

#Write to bams
if(!dir.exists('./filtered_bams/')) dir.create('./filtered_bams/')
selected_reads<-as.logical(rep(0,length(reads[[1]]$seq)));
selected_reads[which(filter_inds)[inds]]<-TRUE;
filterBam(file=paste(bamfdir,bamfname,sep=''),index=baifname,
					destination=paste('./filtered_bams/',sampname,'_filteredreads.bam',sep=''),
					filter=selected_reads,param=params); 
selected_reads<-as.logical(rep(0,length(reads[[1]]$seq)));
selected_reads[which(varread_inds)[inds]]<-TRUE;
filterBam(file=paste(bamfdir,bamfname,sep=''),index=baifname,
					destination=paste('./filtered_bams/',sampname,'_variantreads.bam',sep=''),
					filter=selected_reads,param=params); 
selected_reads<-reads[[1]]$qname%in%var_table$name[
	(!is.na(as.logical(var_table$del_in_target))&as.logical(var_table$del_in_target))|
		(!is.na(as.logical(var_table$ins_in_target))&as.logical(var_table$ins_in_target))];
if(sum(selected_reads)>0){
	filterBam(file=paste(bamfdir,bamfname,sep=''),index=baifname,
						destination=paste('./filtered_bams/',sampname,'_targetvariantreads.bam',sep=''),
						filter=selected_reads,param=params); 
}else{
	print('No variants found in target.')
}
#And all exact matches (useful for some workflows)
selected_reads<-as.logical(rep(0,length(reads[[1]]$seq)));
selected_reads<-reads[[1]]$qname%in%var_table$name[
	(!is.na(as.logical(var_table$del_in_target))&!as.logical(var_table$del_in_target))|
		(!is.na(as.logical(var_table$ins_in_target))&!as.logical(var_table$ins_in_target))];
if(sum(selected_reads)>0){
	filterBam(file=paste(bamfdir,bamfname,sep=''),index=baifname,
						destination=paste('./filtered_bams/',sampname,'_reads_no_tgt_variants.bam',sep=''),
						filter=selected_reads,param=params); 
}else{
	print('No variants found in target.')
}

#Clean up
file.remove(baifname)



