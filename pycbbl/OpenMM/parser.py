def readOpenMMReporterFile(reporter_file):
    with open(reporter_file, 'r') as ef:
        lines = ef.readlines()
        data = {}
        for r in lines[0].split(','):
            data[r.replace('#','').replace('"','').strip()] = []
        for i,r in enumerate(data):
            for line in lines[1:]:
                data[r].append(float(line.strip().split(',')[i]))
    return data
