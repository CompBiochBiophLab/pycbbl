from scrapy.crawler import CrawlerProcess
from .spiders import pdbSpider
import os
import json

def getPDBInformation(pdb_ids, output_file, overwrite=False, download_delay=0.2):
    """
    Retrieves information for a list of PDB ids. The PDB entry page is scrapped
    to extract information that is later stored into a dictionary. The dictionary
    is written into a json file.

    This function uses scrapy as the engine to retrieve information from PDB,
    therefore is not restartable. If the json file already exists, the function
    will load the information from it and return it without executing the scrapy
    spider. This can be changed by giving the option overwrite=True.

    If too many requests are given in an amount of time the site will raise a too
    many request error. The unit of time per request can be controled with the
    download_delay option. Try bigger values if you are experiencing this kind
    error. Otherwise, if you have fewer pdb codes, try smaller values to speed up
    the scrapying.

    Attributes
    ----------
    pdb_ids : list
        List of PDB ids to retrieve information.
    output_file : str
        Path for the output dictionary storing the retrieved information.
    download_delay : float
        Time delay between scrapy requests.
    Returns
    -------
    pdb_data : dict
        Dictionary containing the PDB data.
    """

    if isinstance(pdb_ids, str):
        pdb_ids = [pdb_ids]

    if os.path.exists(output_file) and not overwrite:
        print('Json file %s found.' % output_file)
        print('Reading information from %s file.' % output_file)

    elif not os.path.exists(output_file) or overwrite:
        # create a crawler process with the specified settings
        process = CrawlerProcess(
                  {'USER_AGENT': 'scrapy',
                   'DOWNLOAD_DELAY': download_delay,
                   'LOG_LEVEL': 'ERROR'})

        process.crawl(pdbSpider, pdb_ids=pdb_ids, output_file=output_file)
        process.start()

    output_file = open(output_file)
    pdb_data = json.load(output_file)
    output_file.close()

    return pdb_data
