'''
This file is meant to assist in the downloading of files from a remote server

'''
import os
import ftplib
from multiprocessing import Pool

def _printProgressBar(bytesrecieved, total_size, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (bytesrecieved / float(total_size)))
    filledLength = int(length * bytesrecieved // total_size)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%%' % (prefix, bar, percent), end = '\r')
    # Print New Line on Complete
    if bytesrecieved == total_size:
        print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
        print()

class _Download_tofile:
    def __init__(self, filename, mode):
        self.f = open(filename, mode)
        self.size = 0
    def write(self,s):
        self.size+=len(s)
        self.f.write(s)
    def close(self):
        self.f.close()

class Download_Engine:
    def __init__(self,host,ftp_targets=[],destination=os.getcwd()):
        self.host=host
        #self.ftp=ftplib.FTP(host,user='anonymous',passwd='anonymous@',timeout=20)
        self.ftp_targets=self._clean_host(ftp_targets)
        self.destination=destination

    def _clean_host(self,ftp_targets):
        ftp_list= [ftp_path.split(self.host)[-1] for ftp_path in ftp_targets]
        return(ftp_list)

    def add_target(self,address):
        target=self._clean_host([address])[0]
        self.ftp_targets.append(target)

    def _download(self,address,printprogress=False,in_memory=False):
        path,file=os.path.split(address)
        try:
            ftp=ftplib.FTP(self.host,user='anonymous',passwd='anonymous@',timeout=20)
        except ftplib.all_errors as e:
            raise(e)
        #ftp.login('anonymous','anonymous@')
        try:
            ftp.cwd(path)
        except ftplib.all_errors as e:
            return(e)
        ftp.sendcmd("TYPE i") #swapping out of ASCII mode
        total=int(ftp.size(file))
        f=_Download_tofile(os.path.join(self.destination,file), 'wb')
        if not printprogress:
            print('Downloading: '+file)
        def callback(data):
            f.write(data)
            if printprogress:
                pfix = file
                _printProgressBar(f.size, total_size=total, prefix = pfix, suffix = 'Download Complete', decimals = 1, length = 30, fill = '█')
        try:
            ftp.retrbinary("RETR "+file, callback)
        except ftplib.all_errors as e:
            return(e)
        f.close()
        ftp.close()

    def download(self):
        for target in self.ftp_targets:
            print("Preparing to download "+target+" to "+self.destination)
            self._download(target,printprogress=True)

    def threaded_download(self,num_thread=4):
        pool = Pool(num_thread)
        pool.map(self._download,self.ftp_targets)

    def listdir(self,directory):
        try:
            ftp=ftplib.FTP(self.host,user='anonymous',passwd='anonymous@',timeout=20)
        except ftplib.all_errors as e:
            raise(e)
        try:
            l=ftp.nlst(directory)
        except ftplib.all_errors as e:
            return(e)
        ftp.close()
        return(l)