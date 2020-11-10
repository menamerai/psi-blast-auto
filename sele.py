import numpy as np
import pandas as pd
from Bio import SeqIO
import selenium
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.firefox.options import Options
from soupparser import parse_results
import warnings
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
        self.driver = selenium.webdriver.Firefox(executable_path="./drivers/geckodriver", options=options)

    def get_functional_proteins_in_df(self):
        try: 
            # Split codons
            proteins = self.AA_chain.split("*")
            # Functional proteins generally has > 20 AA (I think)
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

    # I don't really know why did I add this def
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

            # The code below is super ugly.

            # I don't know what happened, but a JavaScript argument[0].submit is not a function
            # error will be raised, and it will stop the python process, but the website still work.
            
            # So I'll have to include this super ugly try except code here,
            # basically ignoring the error.

            try:
                sequence_box.submit()
            except:
                pass
            try:
                WebDriverWait(driver, 60).until(
                    EC.presence_of_element_located((By.CLASS_NAME, "aboveThreshold"))
                )
                result_text = parse_results(driver.current_url)
                result_text.insert(0, protein)
                print(result_text[2] + "\nPositive: " + result_text[4] + "%")
                result_tables.append(result_text)
            except selenium.common.exceptions.TimeoutException:
                # Sometimes, the PSI-BLAST does not work, rasing format errors. 
                # However, weaving two WebDriverWait together to differenciate two results is hard.
                # And the lack of general navigation features I found is quite prominent,
                # so I'll sacrifice some efficiency for the ability to not have a headache for now.
                print("Search timed out, or failed to complete. Skipping amino acid chain " + protein)
                result_tables.append([protein, None, None, None, None])
            count += 1

        driver.quit()
        return pd.DataFrame(result_tables, columns=["Ameno Acid Chain", "DB:ID", "Source", "Identities", "Positives"])