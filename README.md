# Analysis of genotype data from 23andMe DNA testing. 
<p> This R script is based on the <a href="http://www.vincebuffalo.com/blog/2012/03/12/using-bioconductor-to-analyze-your-23andme-data.html">blog post by Vince Buffalo </a>. Some of his codes were outdated, so new ones are used here.The basic workflow involves Bioconductor libraries such as <i>TxDb.Hsapiens.UCSC.hg18.knownGene</i> for transcript annotation, and <i>org.Hs.eg.db</i> for converting annotated Entrez gene IDs to actual gene names, which is helpful for browsing the reference database.</p>
<p>After getting the database of gene names, I queried for the gene symbol "ADH1B" - the protein-coding gene for Alcohol Dehydrogenase II. This enzyme, or the lack of, is the reason behind my flushing when I consume alcohol. From my 23andMe genotype data, I obtaind three hits on Chromosome 4 with the corresponding rsid and genotypes. Next, metadata from the library <i> gwascat</i> was used to merge with my genotype data. This creats an instance with the object "Strongest.SNP.Risk.Allele" showing which mutated base nucleotide is most risky(high correlation with disease in GWAS study). My own risk was calculated by looping over every line of the data table (mapply()) and see whether or not the strongest risk allele mutation is in my genotype. This creates a binary column (False, True) for my risk factors. Since I have hypothroidism, I found the rsid index for anything that matched hypothyroidism under disease.trait in the calculated risk table, and looked at the allele frequency for the transcripts related to hypothyroidism in general. Within the filtered data for hypothyroidsim (13 samples), there wasn't any studies specific to chinese or asian population.</p>
<p>Here're some interesting insights for the "Diease Trait" that showed up in my calculated risk table:</p>
<ul><li>"Economic and political preferences (feminism/equality)"</li>
<li>"Economic and political preferences (environmentalism)"</li>
</ul>
<p> Hmmm, why are they listed under "Disease.Trait". I'm trying to make sense of the fact that DNA testing can predict Economic and political preferences? 😂😂😂</p>

<p> Just for kick, I also plotted a karyogram of my "risk factor" loci to get a broad view of the distribution of strong risk mutations. </p>




<p align="center"><img src='DNA karyograms.jpeg', width=70%, height=70%></p>
<p>note: Just a quick background on 23andMe genetic testing, I received my testing before the FDA shutdown and the regulation changes that followed. Prior to the FDA shut down in November 2013, the company had offered consumers genetic testing to estimate their risk for 240 health conditions (that original report is no longer posted in my account). However, in October 2015 the company received FDA approval to test for carrier genes of 36 diseases that can be passed on to the offsprings.</p>
