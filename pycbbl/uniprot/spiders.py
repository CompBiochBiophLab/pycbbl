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
    allowed_domains = ['www.uniprot.org']

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

        # Save scraped url
        self.uniprot_data[current]['url'] = 'https://www.uniprot.org/uniprot/'+current

        ### Scrape uniprot data here ###

        ## Basic information
        yield scrapy.Request(self.uniprot_data[current]['url'],
                             self.parseBasicInformation,
                             meta={'current': current},
                             dont_filter=True)

        ## Names & Taxonomy
        yield scrapy.Request(self.uniprot_data[current]['url'],
                             self.parseNamesAndTaxonomy,
                             meta={'current': current},
                             dont_filter=True)

        ## Structure ##
        yield scrapy.Request(self.uniprot_data[current]['url'],
                             self.parseStructure,
                             meta={'current': current},
                             dont_filter=True)

        ## Sequence ##
        yield scrapy.Request(self.uniprot_data[current]['url']+'.fasta',
                             self.parseSequence,
                             meta={'current': current},
                             dont_filter=True)


        #
        # ## Funtion information ##
        # self.uniprot_data[current]['Function'] = {}
        #
        # # Sites
        # feature = []
        # for x in response.css('#sitesAnno_section > tr > td > span.context-help::text').getall():
        #     feature.append(x)
        # position = []
        # for x in response.css('#sitesAnno_section > tr > td > a.position::text').getall():
        #     position.append(x)
        # description = []
        # for x in response.css('#sitesAnno_section > tr > td.featdescription > span::text').getall():
        #     description.append(x)
        #
        # if feature != []:
        #     self.uniprot_data[current]['Function']['Sites'] = {}
        #     for i in range(len(feature)):
        #         self.uniprot_data[current]['Function']['Sites'][i+1] = {}
        #         self.uniprot_data[current]['Function']['Sites'][i+1]['Feature key'] = feature[i]
        #         self.uniprot_data[current]['Function']['Sites'][i+1]['Positions'] = position[i]
        #         self.uniprot_data[current]['Function']['Sites'][i+1]['Description'] = description[i]
        #
        # # Go terms
        # self.uniprot_data[current]['Function']['GO - Molecular function'] = []
        # for x in response.css('#function > ul.noNumbering.molecular_function > li > a::attr(href)').getall():
        #     self.uniprot_data[current]['Function']['GO - Molecular function'].append(x.split(':')[-1])
        #
        # self.uniprot_data[current]['Function']['GO - Biological process'] = []
        # for x in response.css('#function > ul.noNumbering.biological_process > li > a::attr(href)').getall():
        #     self.uniprot_data[current]['Function']['GO - Biological process'].append(x.split(':')[-1])
        #
        # ## Family & Domains ##
        #
        # # Domains
        # feature = []
        # for x in response.css('#domainsAnno_section > tr > td > span.context-help::text').getall():
        #     feature.append(x)
        # position = []
        # for x in response.css('#domainsAnno_section > tr > td > a.position::text').getall():
        #     position.append(x.replace(u'\xa0', u' ').replace('–','-'))
        # description = []
        # for x in response.css('#domainsAnno_section > tr > td.featdescription > span::text').getall():
        #     if x != ' ':
        #         description.append(x.replace(u'\xa0', u' '))
        #
        # if feature != []:
        #     self.uniprot_data[current]['Family & Domains'] = {}
        #     if 'Domain' in feature:
        #         self.uniprot_data[current]['Family & Domains']['Domains'] = {}
        #     if 'Repeat' in feature:
        #         self.uniprot_data[current]['Family & Domains']['Repeats'] = {}
        #
        #     dc = 0
        #     rc = 0
        #     for i in range(len(feature)):
        #         if feature[i] == 'Domain':
        #             dc += 1
        #             self.uniprot_data[current]['Family & Domains']['Domains'][dc] = {}
        #             self.uniprot_data[current]['Family & Domains']['Domains'][i+1]['Positions'] = position[i]
        #             self.uniprot_data[current]['Family & Domains']['Domains'][i+1]['Description'] = description[i]
        #         elif feature[i] == 'Repeat':
        #             rc += 1
        #             self.uniprot_data[current]['Family & Domains']['Repeats'][rc] = {}
        #             self.uniprot_data[current]['Family & Domains']['Repeats'][i+1]['Positions'] = position[i]
        #             self.uniprot_data[current]['Family & Domains']['Repeats'][i+1]['Description'] = description[i]
        #

        #


        ### End scraping uniprot data ###

    def parseBasicInformation(self, response):
        current = response.meta['current']

        ## Basic entry information ##
        protein = response.css('#content-protein > h1::text').extract_first()
        self.uniprot_data[current]['Protein'] = protein

        gene = response.css('#content-gene > h2::text').extract_first()
        self.uniprot_data[current]['Gene'] = gene

        organism = response.css('#content-organism > em::text').extract_first()
        self.uniprot_data[current]['Organism'] = organism

    def parseNamesAndTaxonomy(self, response):
        current = response.meta['current']
        self.uniprot_data[current]['Taxonomy'] = []
        #'div#names_and_taxonomy > table > tbody > tr > td'
        for x in response.css('div#names_and_taxonomy > table > tr > td > *'):
            element = x.get().split()[0].split('>')[0].replace('<','')
            if element == 'span':
                tax = x.css('a::text').get()
                if tax != None and tax != 'More...':
                    if 'Ensembl:' not in tax:
                        self.uniprot_data[current]['Taxonomy'].append(tax)
            elif element == 'a':
                if '/taxonomy/' in x.css('a::attr(href)').extract_first():
                    tax = x.css('a::text').extract_first()
                    if not self._isInteger(tax):
                        self.uniprot_data[current]['Taxonomy'].append(tax)

    def parseStructure(self, response):
        current = response.meta['current']
        sequence = ''
        # PDB structures
        self.uniprot_data[current]['Structures'] = []
        for x in response.css('tr > td > a.pdb::text').getall():
            self.uniprot_data[current]['Structures'].append(x)

    def parseSequence(self, response):
        current = response.meta['current']
        sequence = ''
        for x in response.css('p::text').getall():
            for line in x.split('\n'):
                if not line.startswith('>'):
                    sequence += line
        self.uniprot_data[current]['Sequence'] = sequence

    def _isInteger(self, s):
        try:
            int(s)
            return True
        except:
            return False

    def closed(self, spider):
        json.dump(self.uniprot_data, self.output_file)
        self.output_file.close()

class similarProteinSpider(scrapy.Spider):
    """
    Scrapy spider to retrieve the sequence of similar proteins to a target protein.

    Attributes
    ----------
    uniprot_ids : list
        List of uniprot ids to retrieve information.
    output_file : str
        Path for the output dictionary storing the retrieved information.
    """
    allowed_domains = ['www.uniprot.org/']

    def __init__(self, uniprot_id=None, output_file=None, percentage=50, **kwargs):

        self.uniprot_id = uniprot_id
        self.percentage = str(percentage)
        self.similar_proteins = []
        self.output_file = open(output_file, 'w')
        if self.uniprot_id == None:
            raise ValueError('You must give a uniprot ID to retrieve the list of\
            similar proteins.')

    def start_requests(self):
        yield scrapy.Request('https://www.uniprot.org/uniprot/'+self.uniprot_id, self.get_links)

    def get_links(self, response):
        #table-fifty > table > tbody > tr:nth-child(6) > td > small > a
        url = ''
        for x in response.css('#similar_proteins > div > table > tbody > tr > td > small > a::attr(href)').getall():
            if str(self.percentage) in x:
                url = 'https://www.uniprot.org'+x+'.list'
        if url != '':
            yield scrapy.Request(url, self.parse_ids, dont_filter=True)
        else:
            print('Similar protein link not found for code '+self.uniprot_id)

    def parse_ids(self, response):
        for list in response.css('p::text').getall():
            for x in list.split():
                self.similar_proteins.append(x)

    def closed(self, spider):
        json.dump(self.similar_proteins, self.output_file)
        self.output_file.close()
