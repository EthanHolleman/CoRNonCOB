class Log():
    '''
    Methods for writing to log file
    
    '''
    def __init__(self, logfile):
        self.logfile = open(logfile, 'w')
        
    def write_string(self, string):
        self.logfile.write('You typed this string: {}\n'.format(string))
        
    def close_log(self):
        self.logfile.close()
