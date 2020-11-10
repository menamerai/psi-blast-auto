import sele

covid = sele.Genome("MN908947.fna")
covid.psi_blast_search().to_csv("covid.csv")