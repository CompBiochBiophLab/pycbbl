from scrapy.crawler import CrawlerProcess
from .spiders import uniprotSpider, similarProteinSpider
import os
import json

def getInformation(uniprot_ids, output_file, overwrite=False):
    """
    This function retrieves information for a list of uniprot ids. The uniprot
    entry page is scrapped to extract information that is later stored into a dictionary.
    The dictionary is written into a json file.

    This function uses scrapy as the engine to retrieve information from uniprot,
    therefore is not restartable. If the json file already exists, the function
    will load the information from it and return it without executing the scrapy
    spider. This can be changed by giving the option overwrite=True.

    Attributes
    ----------
    uniprot_ids : list
        List of uniprot ids to retrieve information.
    output_file : str
        Path for the output dictionary storing the retrieved information.

    Returns
    -------
    uniprot_data : dict
        Dictionary containing the uniprot data.
    """

    if isinstance(uniprot_ids, str):
        uniprot_ids = [uniprot_ids]

    if os.path.exists(output_file) and not overwrite:
        print('Json file %s found.' % output_file)
        print('Reading information from %s file.' % output_file)

    elif not os.path.exists(output_file) or overwrite:
        # create a crawler process with the specified settings
        process = CrawlerProcess(
                  {'USER_AGENT': 'scrapy',
                   'LOG_LEVEL': 'ERROR'})

        process.crawl(uniprotSpider, uniprot_ids=uniprot_ids, output_file=output_file)
        process.start()

    output_file = open(output_file)
    uniprot_data = json.load(output_file)
    output_file.close()

    return uniprot_data

def getSimilarProteins(uniprot_id, output_file, percentage=50, overwrite=False):
    """
    Gather similar proteins to the target protein appearing in the uniprot database.

    Attributes
    ----------
    uniprot_id : str
        Uniprot ID of the target protein.
    """
    if percentage not in [50, 90]:
        raise ValueError('Only 50 or 90 percentages are valid for uniprot lists of similar porteins.')

    if os.path.exists(output_file) and not overwrite:
        print('Json file %s found.' % output_file)
        print('Reading information from %s file.' % output_file)

    elif not os.path.exists(output_file) or overwrite:
        # create a crawler process with the specified settings
        process = CrawlerProcess(
                  {'USER_AGENT': 'scrapy',
                   'LOG_LEVEL': 'ERROR'})

        process.crawl(similarProteinSpider,
                      uniprot_id=uniprot_id,
                      output_file=output_file, 
                      percentage=percentage)
        process.start()

    output_file = open(output_file)
    uniprot_data = json.load(output_file)
    output_file.close()

    return uniprot_data
