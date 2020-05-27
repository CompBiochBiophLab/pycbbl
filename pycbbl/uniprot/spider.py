import scrapy
import json

class uniprotSpider(scrapy.Spider):
    """
    Scrapy spider to retrieve information for a list of uniprot ids. The uniprot
    entry page is scrapped to extract information that is stored into a dictionary.
    The spider writes this dictionary into a json file.

    Attributes
    ----------
    uniprot_ids : list
        List of uniprot ids to retrieve information.
    output_file : str
        Path for the output dictionary storing the retrieved information.
    """
    allowed_domains = ['www.uniprot.org/']

    def __init__(self, uniprot_ids=None, output_file=None, **kwargs):
        self.uniprot_ids = uniprot_ids
        self.uniprot_data = {}
        self.output_file = open(output_file, 'w')
        if self.uniprot_ids == None:
            raise ValueError('You must give a list with the uniprot IDs to retrieve\
                  information from.')

    def start_requests(self):
        for upid in self.uniprot_ids:
            yield scrapy.Request('https://www.uniprot.org/uniprot/'+upid, self.parse, meta={'upid': upid})

    def parse(self, response):

        # Get input uniprot id
        current = response.meta['upid']
        self.uniprot_data[current] = {}

        ## Scrape uniprot data here ##

        # Basic entry information
        protein = response.css('#content-protein > h1::text').extract_first()
        self.uniprot_data[current]['Protein'] = protein

        gene = response.css('#content-gene > h2::text').extract_first()
        self.uniprot_data[current]['Gene'] = gene

        organism = response.css('#content-organism > em::text').extract_first()
        self.uniprot_data[current]['Organism'] = organism

        # Go terms
        self.uniprot_data[current]['GO - Molecular function'] = []
        for x in response.css('#function > ul.noNumbering.molecular_function > li > a::attr(href)').getall():
            self.uniprot_data[current]['GO - Molecular function'].append(x.split(':')[-1])

        self.uniprot_data[current]['GO - Biological process'] = []
        for x in response.css('#function > ul.noNumbering.biological_process > li > a::attr(href)').getall():
            self.uniprot_data[current]['GO - Biological process'].append(x.split(':')[-1])

        # PDB structures
        self.uniprot_data[current]['Structures'] = []
        for x in response.css('tr > td > a.pdb::text').getall():
            self.uniprot_data[current]['Structures'].append(x)

        ## End scraping uniprot data ##

    def closed(self, spider):
        json.dump(self.uniprot_data, self.output_file)
        self.output_file.close()
