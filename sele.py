import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from selenium import webdriver
import selenium
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.firefox.options import Options
import warnings
import soupparser
warnings.simplefilter(action='ignore')


class Genome:

    def __init__(self, filepath):
        self.filepath = filepath
        self.DNA = SeqIO.read(filepath, "fasta").seq
        self.mRNA = self.DNA.transcribe()
        self.AA_chain = self.mRNA.translate()
        self.func_df = None

        options = Options()
        options.headless = True
        self.driver = webdriver.Firefox(executable_path="./drivers/geckodriver", options=options)

    def get_functional_proteins_in_df(self):
        try: 
            proteins = self.AA_chain.split("*")
            func_proteins = list(
                filter(lambda protein: len(protein) > 20, proteins))
            func_df = pd.DataFrame({"Functional Proteins": func_proteins})
            func_df["Functional Proteins"] = func_df["Functional Proteins"].apply(
                lambda x: "".join(x))
            if self.func_df == None:
                self.func_df = func_df
            return func_df
        except:
            print("An error occured. (get_functional_proteins)")

    def update_functional_proteins(self, df=None):
        if not isinstance(df, pd.DataFrame):
            raise ValueError("The argument passes must be a dataframe!")
        try: 
            self.func_df = df
            print("Updated.")
            print(self.func_df.head())
        except:
            print("An error occured. Please check if the dataframe is valid.")
            print("Attempting to update dataframe via class data instead...")
            self.func_df = self.get_functional_proteins_in_df()

    
    def psi_blast_search(self):
        if self.func_df == None:
            print("A functional proteins dataframe has not been created.\nAttempting to create one...")
            self.func_df = self.get_functional_proteins_in_df()
        else:
            print("Existing dataframe detected.")
        driver = self.driver
        protein_list = self.func_df["Functional Proteins"].values.tolist()
        result_tables = []
        count = 1 

        for protein in protein_list:
            print("=" * 10)
            print(str(count) + ". " + protein)
            driver.get("https://www.ebi.ac.uk/Tools/sss/psiblast/")
            sequence_box = driver.find_element_by_id("sequence")
            sequence_box.send_keys(protein)
            try:
                sequence_box.submit()
            except:
                pass
            try:
                WebDriverWait(driver, 60).until(
                    EC.presence_of_element_located((By.CLASS_NAME, "aboveThreshold"))
                )
                result_text = soupparser.parse_results(driver.current_url)
                result_text.insert(0, protein)
                print(result_text[2])
                result_tables.append(result_text)
            except selenium.common.exceptions.TimeoutException:
                print("Search timed out, or failed to complete. Skipping amino acid chain " + protein)
                result_tables.append([protein, None, None, None, None])
            count += 1

        driver.quit()
        return pd.DataFrame(result_tables, columns=["Ameno Acid Chain", "DB:ID", "Source", "Identities", "Positives"])