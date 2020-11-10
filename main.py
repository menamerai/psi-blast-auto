from sele import Genome

covid = Genome("MN908947.fna")
covid.psi_blast_search().to_csv("covid.csv")