from bs4 import BeautifulSoup
import requests

def parse_results(url):
    page = requests.get(url)
    soup = BeautifulSoup(page.content, "html.parser")

    results = soup.find_all("td")
    for i in range(len(results)):
        results[i] = results[i].get_text().rstrip("\n").strip()

    return [results[1], results[2], results[5], results[6]]
