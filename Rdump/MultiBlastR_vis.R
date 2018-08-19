##########
#	This is a script to visualise the output drom blastfest
#
# For single fasta file
# sh /Users/sa01fd/Github/e-utils/MultiBlastP.sh -f /Users/sa01fd/testUltra/glycolisis_kegg_salarius.fasta 
#
# running on many fasta files in batch
# for f in ./fastas/*.fasta; do sh /Users/sa01fd/Github/e-utils/MultiBlastP.sh -f $f; done
##########

RAxMLRinke = read.tree("~/pruned_tree")
RAxMLANVIORooted =  midpoint.root(RAxMLRinke)
tree2 <- ladderize(RAxMLANVIORooted, right = FALSE)

is_tip <- tree2$edge[,2] <= length(tree2$tip.label)
ordered_tips <- tree2$edge[is_tip, 2]
tree2$tip.label[ordered_tips]


selected='Phenylacetic_acid'
selected='benzyl_alcohol'
selected='kegg_denitri_VT8'
selected='glycolisis_kegg_salarius'
selected='4-hydroxy_benzoate'
selected='4-Hydroxy_phenylacetate'
selected='Phenyl_propanoid_'
selected='transport_psychro'
selected='sulfur_psychro'
selected='nitrogen_VT8'
selected='chemotaxis_HP15'
selected='PTS'
selected='some_vits'
selected='glycogen_metabolism'
selected='D.shibae_cobalamin'

selected='DG893_bacteriocin_antismash'
selected='DshiDMSO'
selected='DshiPhoto'
Dshi_NADH_mnhOperon
DshiDMSO
selected='Ethanolamine_MAQU'


selected='HP15PTS1'
selected='HP15PTS2'

selected='NKSG1_ATPase'
selected='NKSG1_ATPase2'
selected='alginate_SK2'
selected='alginate_pseudo'
selected='catechol_meta'

selected='catechol_ortho'
selected='cbb3'
selected='cobalamin'
selected='cobalaminATCC'
selected='cobalaminATCC'

selected='coenzsmeF390_HP15'
selected='cytoC'
extracellular_poly
selected='extracellular_poly'
2nuo_pseudomonas
selected='2nuo_pseudomonas'


selected='formate1_dehydro_NKSG1'
selected='formate_pseudomonas'
selected='cobalaminATCC'
full_glycolate
selected='Conjugal_plasmid2'

selected='glucose_mannose'
selected='heme_exporter_HP15'
selected='heme_exporter_HP15'
selected='lactate_cluster'
selected='lipid_A_cluster'
selected='maltose'
selected='maltose_redo'
selected='multi_subunit'
selected='murein_cluster2'
selected='ndh2_pseudo'
selected='phn_operon_psuedo'

selected='extracellular_poly'
selected='TRAP_mannitol'
selected='betacarotene_test'
selected='PTS_DG893'


df = read.table(paste('~/testUltra/',selected,'/',selected,'.abundance_table.txt',sep=''),header=TRUE)
df.annot = read.table(paste('~/testUltra/',selected,'/',selected,'.annotation_file.txt',sep=''),sep='\t')


rownames(df.annot)=df.annot$V1

df = data.frame(df)
df[is.na(df)] <- 0
df[df>1] <- 1

#this is needed to tidy up the genome names and make them coherent (if you were smarter, you can ignore this)
colnames(df)=gsub("\\.","_",colnames(df))
colnames(df)=gsub("_SW_145","",colnames(df))
colnames(df)=gsub("Marinobacter_salexigens_HJR7","Marinobacter_salexigens_strain_HJR7",colnames(df))
colnames(df)=gsub("Marinobacter_guineae_M3B","Marinobacter_guineae_strain_M3B",colnames(df))
colnames(df)=gsub("Marinobacter_adhaerens_PBVC038","Marinobacter_adhaerens_strain_PBVC038",colnames(df))




df = df[,intersect(tree2$tip.label[ordered_tips],colnames(df))]
df$gene = rownames(df)
df$gene2 = df.annot[rownames(df),'V2']

df1 = melt(df,id=c('gene','gene2'))

#levels(df1$variable) =  tree2$tip.label[ordered_tips]

#df1 $variable <- factor(df1 $variable, levels = RAxMLANVIORooted$tip.label)


ggplot(df1, aes(y = variable,x = gene2)) + 
     geom_point(aes(fill = value),colour="black", size =2.5, shape=21)  +theme_bw() +
     labs(x=NULL, y = NULL)+scale_fill_gradient2(low = "white", high = "black", midpoint = 0)+ggtitle(selected)+theme(axis.text=element_text(size=6),
        axis.title=element_text(size=6))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        ggplot(df1, aes(y = variable,x = gene2)) + 
     geom_point(aes(fill = value),colour="black", size =2.5, shape=21)  +theme_bw() +
     labs(x=NULL, y = NULL)+scale_fill_gradient2(low = "white", high = "black", midpoint = 0)+ggtitle(selected)+theme(axis.text=element_text(size=1),
        axis.title=element_text(size=1))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
       