import scrapy
import json
import time

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
        self.parseBasicInformation(response, current)

        ## Function
        self.parseFunction(response, current)

        ## Names & Taxonomy
        self.parseNamesAndTaxonomy(response, current)

        ## Structure ##
        self.parseStructure(response, current)

        ## GO terms ##
        self.parseGO(response, current)

        ## Sequence ##
        yield scrapy.Request(self.uniprot_data[current]['url']+'.fasta',
                             self.parseSequence,
                             meta={'current': current},
                             dont_filter=True)

        ### Finish scraping uniprot data here ###

    def parseFunction(self, response, current):
        self.uniprot_data[current]['Function'] = {}

        function = response.css('div#function')

        # Catalysis

        kcat = function.xpath('//span[contains(@title, "Michaelis")]/../../../span/text()').extract()
        Kms = []
        for Km in function.xpath('//span[contains(@title, "Michaelis")]/../span[last()]/text()').extract():
            Kms.append(Km)

        if Kms != []:
            self.uniprot_data[current]['Function']['Catalysis'] = {}
            self.uniprot_data[current]['Function']['Catalysis']['Km'] = Kms
        if kcat != []:
            self.uniprot_data[current]['Function']['Catalysis']['Kcat'] = kcat[0]

        # Sites
        feature_keys = []
        positions = []
        descriptions = []

        for i,x in enumerate(function.css('#sitesAnno_section td > :nth-child(1)::text').extract()):

            if i%3 == 0:
                feature_keys.append(x)
            elif i%3 == 1:
                positions.append(x)
            elif i%3 == 2:
                descriptions.append(x)

        self.uniprot_data[current]['Function']['Sites'] = {}
        sites = self.uniprot_data[current]['Function']['Sites']

        for i,x in enumerate(positions):
            sites[x] = {}
            sites[x]['Feature key'] = feature_keys[i]
            sites[x]['Description'] = descriptions[i]

    def parseBasicInformation(self, response, current):

        ## Basic entry information ##
        protein = response.css('#content-protein > h1::text').extract_first()
        self.uniprot_data[current]['Protein'] = protein

        gene = response.css('#content-gene > h2::text').extract_first()
        self.uniprot_data[current]['Gene'] = gene

        organism = response.css('#content-organism > em::text').extract_first()
        self.uniprot_data[current]['Organism'] = organism

    def parseNamesAndTaxonomy(self, response, current):

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

    def parseStructure(self, response, current):
        sequence = ''
        # PDB structures
        self.uniprot_data[current]['Structures'] = []
        for x in response.css('tr > td > a.pdb::text').getall():
            self.uniprot_data[current]['Structures'].append(x)

    def parseGO(self, response, current):
        self.uniprot_data[current]['GO'] = {}

        self.uniprot_data[current]['GO']['Molecular Function'] = []
        for x in response.css('#function .molecular_function > li > a::attr(href)').getall():
            self.uniprot_data[current]['GO']['Molecular Function'].append(x.split('GO:')[-1])

        self.uniprot_data[current]['GO']['Biological Process'] = []
        for x in response.css('#function .biological_process > li > a::attr(href)').getall():
            self.uniprot_data[current]['GO']['Biological Process'].append(x.split('GO:')[-1])

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

    def __init__(self, uniprot_ids=None, output_file=None, percentage=50, **kwargs):

        if isinstance(uniprot_ids, str):
            self.uniprot_ids = [uniprot_ids]
        elif isinstance(uniprot_ids, list):
            self.uniprot_ids = uniprot_ids
        self.percentage = str(percentage)
        self.similar_proteins = {}
        self.output_file = open(output_file, 'w')
        if self.uniprot_ids == None:
            raise ValueError('You must give a list of uniprot IDs to retrieve the list of\
            similar proteins.')

    def start_requests(self):
        for upid in self.uniprot_ids:
            yield scrapy.Request('https://www.uniprot.org/uniprot/'+upid, self.get_links, meta={'upid': upid})

    def get_links(self, response):
        #table-fifty > table > tbody > tr:nth-child(6) > td > small > a
        url = ''
        for x in response.css('#similar_proteins > div > table > tbody > tr > td > small > a::attr(href)').getall():
            if str(self.percentage) in x:
                url = 'https://www.uniprot.org'+x+'.list'
        if url != '':
            yield scrapy.Request(url, self.parse_ids, meta={'upid': response.meta['upid']}, dont_filter=True)
        else:
            print('Similar protein link not found for code '+self.uniprot_id)

    def parse_ids(self, response):
        current = response.meta['upid']
        self.similar_proteins[current] = {}
        self.similar_proteins[current]['codes'] = []
        self.similar_proteins[current]['sequences'] = {}

        for list in response.css('p::text').extract():
            for x in list.split():
                self.similar_proteins[current]['codes'].append(x)

        # Retrieve sequences of similar proteins
        for s in self.similar_proteins[current]['codes']:
            url = 'https://www.uniprot.org/uniprot/'+s+'.fasta'
            yield scrapy.Request(url, self.getSequences,
                                 meta={'current': current,
                                       'upid': s},
                                       dont_filter=True)

    def getSequences(self, response):
        current = response.meta['current']
        upid = response.meta['upid']
        sequence = ''
        for x in response.css('p::text').getall():
            for line in x.split('\n'):
                if not line.startswith('>'):
                    sequence += line
        self.similar_proteins[current]['sequences'][upid] = sequence

    def closed(self, spider):
        json.dump(self.similar_proteins, self.output_file)
        self.output_file.close()

class pdbSpider(scrapy.Spider):
    """
    Scrapy spider to retrieve information for a list of PDB ids. The PDB
    entry page is scrapped to extract information that is stored into a dictionary.
    The spider writes this dictionary into a json file.

    Attributes
    ----------
    pdb_ids : list
        List of uniprot ids to retrieve information.
    output_file : str
        Path for the output dictionary storing the retrieved information.
    """
    allowed_domains = ['www.rcsb.org/']

    def __init__(self, pdb_ids=None, output_file=None, **kwargs):
        self.pdb_ids = pdb_ids
        self.pdb_data = {}
        self.output_file = open(output_file, 'w')
        if self.pdb_ids == None:
            raise ValueError('You must give a list with the PDB IDs to retrieve\
                  information.')

    def start_requests(self):
        for pdbid in self.pdb_ids:
            yield scrapy.Request('https://www.rcsb.org/structure/'+pdbid, self.parse, meta={'pdbid': pdbid})

    def parse(self, response):

        # Get input uniprot id
        current = response.meta['pdbid']
        self.pdb_data[current] = {}

        # Save scraped url
        self.pdb_data[current]['url'] = 'https://www.rcsb.org/structure/'+current

        ### Scrape PDB data here ###

        ## Basic information
        self.parseBasicInformation(response, current)

        ## Macromolecules information
        self.parseMacromolecules(response, current)

        ## Small Molecules information
        self.parseSmallMolecules(response, current)

    def parseBasicInformation(self, response, current):
        ## Basic entry information ##
        structureTitle = response.css('#structureTitle::text').extract_first()
        self.pdb_data[current]['Title'] = structureTitle

        # Experimental data
        # for x in css.response('#exp_header_0_snapshot li').getall():
            # print(x)

        resolution = response.css('#exp_header_0_diffraction_resolution::text').extract_first()
        if resolution != None:
            resolution = float(resolution.replace('Å',''))
        self.pdb_data[current]['Resolution'] = resolution

    def parseMacromolecules(self, response, current):
        ## Macromolecules entry information ##
        self.pdb_data[current]['Macromolecules'] = {}
        molecule_names = []
        chains = []
        organisms = []
        lengths = []
        uniprot_ids = []

        # Scrape macromolecules table
        for i,x in enumerate(response.css('#MacromoleculeTable tr:nth-child(3) td')):
            index = (i+1)%6
            if index== 1:
                molecule_name = x.css('td::text').extract_first()
                molecule_names.append(molecule_name)
            elif index == 2:
                chain_ids = x.css('a::text').extract()
                chains.append(chain_ids)
            elif index == 3:
                length = x.css('td::text').extract_first()
                lengths.append(length)
            elif index == 4:
                organism = x.css('td a::text').extract_first()
                organisms.append(organism)

        for i,x in enumerate(response.css('table')):
            id_data = x.css('::attr(id)').extract_first()
            if id_data != None and 'table_macromolecule-protein-entityId' in id_data:
                upid = x.css('.text-left .querySearchLink::text').extract()
                if upid == []:
                    uniprot_ids.append(None)
                else:
                    uniprot_ids.append(upid[0])

        # Store data in dictionary
        for entity in zip(chains, molecule_names, lengths, organisms, uniprot_ids):
            for chain in entity[0]:
                self.pdb_data[current]['Macromolecules'][chain] = {}
                self.pdb_data[current]['Macromolecules'][chain]['Molecule'] = entity[1]
                self.pdb_data[current]['Macromolecules'][chain]['Length'] = entity[2]
                self.pdb_data[current]['Macromolecules'][chain]['Organism'] = entity[3]
                self.pdb_data[current]['Macromolecules'][chain]['Uniprot'] = entity[4]

    def parseSmallMolecules(self, response, current):
        ## Small Molecules entry information ##
        self.pdb_data[current]['Small Molecules'] = {}
        ligand_ids = []
        chains = []
        names = []
        for i,x in enumerate(response.css('#LigandsMainTable td')):
            if (i+1)%5 == 1:
                ligand_ids.append(x.css('a::text').extract_first())
            elif (i+1)%5 == 2:
                chains.append(x.css('div.ellipsisToolTip::text').extract())
            elif (i+1)%5 == 3:
                names.append(x.css('strong::text').extract_first())

        for chain in chains:
            for c in chain:
                self.pdb_data[current]['Small Molecules'][c] = {}
                self.pdb_data[current]['Small Molecules'][c]['ID'] = []
                self.pdb_data[current]['Small Molecules'][c]['Name'] = []

        for entity in zip(chains, ligand_ids, names):
            for chain in entity[0]:
                self.pdb_data[current]['Small Molecules'][chain]['ID'].append(entity[1])
                self.pdb_data[current]['Small Molecules'][chain]['Name'].append(entity[2])

    def closed(self, spider):
        json.dump(self.pdb_data, self.output_file)
        self.output_file.close()
