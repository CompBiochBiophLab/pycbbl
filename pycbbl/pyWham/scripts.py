import xml.etree.ElementTree as ET
from xml.dom import minidom

class pyWhamXmlFile:

    def __init__(self):

        self.xml = ET
        self.root = self.xml.Element('WhamSpec')
        self.general = self.xml.SubElement(self.root, 'General')
        self.coordinates = self.xml.SubElement(self.general, 'Coordinates')
        self.defaultCoordinateFileReader = self.xml.SubElement(self.general, 'DefaultCoordinateFileReader')
        self.binnings = self.xml.SubElement(self.general, 'Binnings')
        self.trajectories = self.xml.SubElement(self.root, 'Trajectories')
        self.jobs = self.xml.SubElement(self.root, 'Jobs')
        self.coordinate = {}
        self.binning = {}
        self.trajectory = {}

    def add_coordinate(self, coordinate_name):

        self.coordinate[coordinate_name] = self.xml.SubElement(self.coordinates, 'Coordinate', name=coordinate_name)

    def add_defaultCoordinateFileReader(self, returns_time=False):
        self.returns_time = returns_time
        self.defaultCoordinateFileReader.set('returnsTime',str(returns_time).lower())
        if returns_time == True:
            self.returnTimeList = self.xml.SubElement(self.defaultCoordinateFileReader, 'ReturnTimeList')
        self.return_list = {}
        for c in self.coordinate:
            self.return_list[c] = self.xml.SubElement(self.defaultCoordinateFileReader, 'ReturnList', name=c)

    def add_binning(self, coordinate_name, interval=None, begin=None, end=None):

        if interval == None:
            raise ValueError('interval = None. You must specify the binning intervalof the coordinate')

        if begin != None and end == None or begin == None and end != None:
            raise ValueError('You must give both being and end statements or None')

        self.binning[coordinate_name] = self.xml.SubElement(self.binnings, 'Binning', name=coordinate_name)
        if begin != None:
            self.binning['begin'] = self.xml.SubElement(self.binning[coordinate_name], 'Begin')
            self.binning['begin'].text = str(begin)
        if end != None:
            self.binning['end'] = self.xml.SubElement(self.binning[coordinate_name], 'End')
            self.binning['end'].text = str(end)
        self.binning['interval'] = self.xml.SubElement(self.binning[coordinate_name], 'Interval')
        self.binning['interval'].text = str(interval)

    def add_trajectory(self, trajectory_file, temperature=None, tBegin=None, energy_function=None):

        if temperature == None:
            raise ValueError('T = None, you need to give a temperature for the given trajectory')
        if energy_function == None:
            raise ValueError('energy_function = None, you need to specify the variable corresponding to the energy of the system')
        if self.returns_time == True and tBegin == None:
            raise ValueError('returns_time = True, but not tBegin was given for the trajectory. Please set tBegin as to when to start reading frames in the given trajectory')

        self.trajectory[trajectory_file] = self.xml.SubElement(self.trajectories, 'Trajectory')
        self.trajectory[trajectory_file].set('T', str(temperature))
        if self.returns_time == True:
            self.trajectory[trajectory_file].set('tBegin', str(tBegin))
        self.energyFunction = self.xml.SubElement(self.trajectory[trajectory_file], 'EnergyFunction')
        self.energyFunction.text = energy_function
        self.coordinateFiles = self.xml.SubElement(self.trajectory[trajectory_file], 'CoordinateFiles')
        self.coordinateFile = self.xml.SubElement(self.coordinateFiles, 'CoordinateFile')
        self.coordinateFile.text = trajectory_file

    # Add Jobs
    def DoS(self, output_file):
        self.dos = self.xml.SubElement(self.jobs, 'DensityOfStates')
        self.dos.set('outFile', output_file)

    def freeEnergy(self, output_file, coordinates=None, energy_function=None, temperature=None, parameters=None):
        if coordinates == None:
            raise ValueError('coordinates = None, you need to give coordinates upon which the free energy will be projected')
        if energy_function == None:
            raise ValueError('energy_function = None, you need to give and energy function to calculate the free energy')
        if temperature == None:
            raise ValueError('temperature = None, you need to give a temperature to calculate the free energy')

        self.freeEnergy = self.xml.SubElement(self.jobs, 'FreeEnergy')
        self.freeEnergy.set('outFilePrefix', output_file)

        self.coordinatesFE = self.xml.SubElement(self.freeEnergy, 'Coordinates')
        if isinstance(coordinates,str):
            self.coordinateFE = self.xml.SubElement(self.coordinatesFE, 'Coordinate')
            self.coordinateFE.set('name', coordinates)
        if isinstance(coordinates,tuple):
            for c in coordinates:
                self.coordinateFE = self.xml.SubElement(self.coordinatesFE, 'Coordinate')
                self.coordinateFE.set('name', c)

        self.energyFunctionFE = self.xml.SubElement(self.freeEnergy, 'EnergyFunction')
        self.energyFunctionFE.text = energy_function

        self.temeperaturesFE = self.xml.SubElement(self.freeEnergy, 'Temperatures')
        self.temeperaturesFE.text = str(temperature)

        if parameters != None:
            self.parametersFE = self.xml.SubElement(self.freeEnergy, 'Parameters')
            if isinstance(parameters, str):
                self.parameterFE = self.xml.SubElement(self.parametersFE, 'Parameter')
                self.parameterFE.set('name', parameters)
                self.parameterFE.text = 'true'
            if isinstance(parameters, tuple) or isinstance(parameters, list):
                for p in parameters:
                    self.parameterFE = self.xml.SubElement(self.parametersFE, 'Parameter')
                    self.parameterFE.set('name', p)
                    self.parameterFE.text = 'true'

    def heatCapacity(self, output_file, energy_function=None, temperatures=None):

        if energy_function == None:
            raise ValueError('energy_function = None, you need to give and energy function to calculate the heat capacity')
        if temperatures == None:
            raise ValueError('temperatures = None, you need to give the temperatures to calculate the heat capacity')

        self.heatCapacity = self.xml.SubElement(self.jobs, 'HeatCapacity')
        self.heatCapacity.set('outFile', output_file)

        self.energyFunctionHC = self.xml.SubElement(self.heatCapacity, 'EnergyFunction')
        self.energyFunctionHC.text = energy_function

        self.temperaturesHC = self.xml.SubElement(self.heatCapacity, 'Temperatures')
        if isinstance(temperatures, str):
            self.temperaturesHC.text = temperatures
        elif isinstance(temperatures, list):
            temperatures = ','.join(temperatures)
            self.temperaturesHC.text = temperatures

    def write_xml(self, file_name):
        xmlstr = minidom.parseString(self.xml.tostring(self.root)).toprettyxml(indent="  ")
        with open(file_name, "w") as f:
            f.write(xmlstr)
