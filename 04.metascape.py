### cmd: java -Dwebdriver.gecko.driver="geckodriver.exe" -jar selenium-server-4.0.0.jar standalone

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.options import Options as FirefoxOptions
import time
import pandas as pd
from urllib.request import urlretrieve
from threading import Thread, Semaphore
import os

def process_column(clusterIdx, geneLst, outfolder, semaphore):
    with semaphore:  # This will block if there are already 8 threads using the semaphore
        options = FirefoxOptions()
        driver = webdriver.Firefox(options=options)
        try:
            driver.get("https://metascape.org/gp/index.html#/main/step1")
            time.sleep(5)
            textArea = driver.find_element("xpath", '/html/body/div[3]/div/div[1]/div[1]/div/table/tbody/tr[1]/td[2]/table/tbody/tr[2]/td/div/textarea')
            textArea.send_keys(geneLst)
            submit = driver.find_element("xpath", '/html/body/div[3]/div/div[1]/div[1]/div/table/tbody/tr[1]/td[2]/table/tbody/tr[4]/td/div/button[1]')
            submit.click()
            analysis = driver.find_element("xpath", '//*[@id="button-div"]/button[1]')
            analysis.click()
            time.sleep(10)
            process = driver.find_element("xpath", '/html/body/div[3]/div/div[1]/div[1]/div/table/tbody/tr[3]/td[2]/div[2]/div[4]/span/span').text
            while float(process) < 100:
                print(process)
                time.sleep(60)
                process = driver.find_element("xpath", '/html/body/div[3]/div/div[1]/div[1]/div/table/tbody/tr[3]/td[2]/div[2]/div[4]/span/span').text
            report_link = driver.find_element("xpath", '/html/body/div[3]/div/div[1]/div[1]/div/table/tbody/tr[3]/td[2]/div[2]/div[3]/a').get_attribute("href")
            sessionID = report_link.split("session_id=")[1].split("&")[0]
            resSite = f"https://metascape.org/gp/index.html#/reportfinal/{sessionID}"
            driver.get(resSite)
            time.sleep(5)
            download_link = driver.find_element("xpath", '/html/body/div[3]/div/div/div/div[2]/a[3]').get_attribute("href")
            outfile = os.path.join(outfolder, f"{clusterIdx}.metascape.zip")
            urlretrieve(download_link, outfile)
        finally:
            driver.quit()

file = "C:/Users/qqiu/Personal Data/Work/Projects/Multiomics/data/strain_wise.DEG.merged.wide.out" # "C:/Users/qqiu/Personal Data/Work/Projects/Multiomics/data/DEG.merged.wide.out"
outfolder = "C:/Users/qqiu/Personal Data/Work/Projects/Multiomics/GO_metascape/"
dat = pd.read_table(file, header=0)

# Create a semaphore that allows up to 8 concurrent threads
semaphore = Semaphore(8)

threads = []
for i in range(1, len(dat.columns)):
    clusterIdx = dat.columns[i].replace("/", ".")
    valid_data = dat[dat.iloc[:, i].notna()].iloc[:, i]
    if len(valid_data) > 30:
        geneLst = "\n".join(valid_data[:3000] if len(valid_data) > 3000 else valid_data)
        thread = Thread(target=process_column, args=(clusterIdx, geneLst, outfolder, semaphore))
        threads.append(thread)
        thread.start()

for thread in threads:
    thread.join()
