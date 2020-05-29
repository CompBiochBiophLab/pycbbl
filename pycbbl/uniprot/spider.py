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

        # Funtion information
        self.uniprot_data[current]['Function'] = {}

        # Sites
        feature = []
        for x in response.css('#sitesAnno_section > tr > td > span.context-help::text').getall():
            #sitesAnno_section > tbody > tr:nth-child(2) > td:nth-child(1) > span
            feature.append(x)
        position = []
        for x in response.css('#sitesAnno_section > tr > td > a.position::text').getall():
            #sitesAnno_section > tbody > tr:nth-child(2) > td:nth-child(1) > span
            position.append(x)
        description = []
        for x in response.css('#sitesAnno_section > tr > td.featdescription > span::text').getall():
            #sitesAnno_section > tbody > tr:nth-child(2) > td:nth-child(1) > span
            description.append(x)

        if feature != []:
            self.uniprot_data[current]['Function']['Sites'] = {}
            for i in range(len(feature)):
                self.uniprot_data[current]['Function']['Sites'][i+1] = {}
                self.uniprot_data[current]['Function']['Sites'][i+1]['Feature key'] = feature[i]
                self.uniprot_data[current]['Function']['Sites'][i+1]['Positions'] = position[i]
                self.uniprot_data[current]['Function']['Sites'][i+1]['Description'] = description[i]

        # Go terms
        self.uniprot_data[current]['Function']['GO - Molecular function'] = []
        for x in response.css('#function > ul.noNumbering.molecular_function > li > a::attr(href)').getall():
            self.uniprot_data[current]['Function']['GO - Molecular function'].append(x.split(':')[-1])

        self.uniprot_data[current]['Function']['GO - Biological process'] = []
        for x in response.css('#function > ul.noNumbering.biological_process > li > a::attr(href)').getall():
            self.uniprot_data[current]['Function']['GO - Biological process'].append(x.split(':')[-1])

        # Sequence
        sequences = []
        sequence = ''
        seq_count = 0
        prev_seq_count = 0
        index = 0
        for x in response.css('pre.sequence::text').getall():
            try:
                seq_count = max([int(c) for c in x.split()])
            except:
                for p in x:
                    if p != ' ':
                        sequence += p
            if seq_count < prev_seq_count:
                sequences.append(sequence)
                sequence = ''
            prev_seq_count = seq_count
        sequences.append(sequence)
        if len(sequences) == 1:
            self.uniprot_data[current]['Sequence'] = sequences[0]
        else:
            for i,seq in enumerate(sequences):
                self.uniprot_data[current]['Sequence'+str(i+1)] = seq

        # PDB structures
        self.uniprot_data[current]['Structures'] = []
        for x in response.css('tr > td > a.pdb::text').getall():
            self.uniprot_data[current]['Structures'].append(x)

        ## End scraping uniprot data ##

    def closed(self, spider):
        json.dump(self.uniprot_data, self.output_file)
        self.output_file.close()
