'''
This file is meant to assist in the downloading of files from a remote server

'''



from threading import Thread
from queue import Queue
import ftplib,gzip,os
from datetime import datetime
class ProgressBar:
    def __init__(self,total_size, prefix = '', suffix = '', decimals = 1, length = 25, fill = '█'):
        self.bytesrecieved=0
        self.total_size = total_size
        self.prefix = prefix
        self.suffix = suffix
        self.decimals = decimals
        self.length = length
        self.fill = fill
    def print(self):
        percent = ("{0:." + str(self.decimals) + "f}").format(100 * (self.bytesrecieved / float(self.total_size)))
        filledLength = int(self.length * self.bytesrecieved // self.total_size)
        bar = self.fill * filledLength + '-' * (self.length - filledLength)
        print('\r%s |%s| %s%%' % (self.prefix, bar, percent), end = '\r')
        # Print New Line on Complete
        if self.bytesrecieved == self.total_size:
            print('\r%s |%s| %s%% %s' % (self.prefix, bar, percent, self.suffix), end = '\r')
            print()
    def recievebytes(self,b):
        self.bytesrecieved+=len(b)

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

class logger:
    def __init__(self):
        self.failed = []
        self.passed =[]
    def append(self,signal_status):
        identity,lcl_file,e = signal_status
        if isinstance(e,str):
            if '226' in e:
                self.passed.append((identity,lcl_file))
            else:self.failed.append((identity,e))
        else:
            self.failed.append((identity,e))
    def __str__(self):
        msg=["Download Engine Logger Report","\tSucessfuly Downloaded"]
        for el in self.passed:
            msg.append('\t\t'+' : '.join(el))
        msg.append("\tDownload Failures and Reason")
        for el in self.failed:
            msg.append('\t\t'+' : '.join(el))
        return('\n'.join(msg))
    def __repr__(self):
        msg="Download_Engine Logger:\n\tComplete = {}\n\tFailed = {}".format(str(len(self.passed)),str(len(self.failed)))
        return (msg)

    def get_failed(self):
        return [el[0] for el in self.failed]

    def get_passed(self):
        return [el[0] for el in self.passed]

    def extend(self,log):
        self.failed.extend(log.failed)
        self.passed.extend(log.passed)
        
    def save_log(self,location):
        fname = os.path.join(location,'Download_{}.log'.format(datetime.now().strftime("%Y-%m-%d_%H:%M")))
        with open(fname,'w') as f:
            for entry in self.passed:
                line = ['passed',entry[0],entry[1]]
                f.write(','.join(line)+'\n')
            for entry in self.failed:
                line = ['failed',entry[0],str(entry[1])]
                f.write(','.join(line)+'\n')
                


class Download_Engine:
    def __init__(self,remote_files=[],destination=None):
        #self.ftp=ftplib.FTP(host,user='anonymous',passwd='anonymous@',timeout=20)

        #URL queue should contain a tuple of (id, host, address)
        self.URL_queue = Queue()
        self._process_remote(remote_files)
        self.logger = logger()

        if destination == None:
            self.destination=os.getcwd()
        else:
            if destination[0] == '~':
                self.destination = os.path.join(os.path.expanduser('~'),destination[2:])
            else:
                self.destination=destination
            if not os.path.exists(self.destination):
                raise FileNotFoundError('The dirctory "{}" does not exsit'.format(self.destination))

    def add_remotefile(self,entry):
        identity,file = entry
        address=file.split('//')[-1].split('/')
        host = address.pop(0)
        self.URL_queue.put((identity,host,'/'.join(address)))

    def _process_remote(self,remote_files):
        for entry in remote_files:
            self.add_remotefile(entry)


    def _touch_files(self):
        while not self.URL_queue.empty():
            report=[]
            refill=[]
            host, file = self.URL_queue.get()
            refill.append((host, file))
            ftp = self._connect(host)
            try:
                lst=ftp.nlst(file)
                if lst:
                    report.append(True)
            except ftplib.error_temp:
                report.append(False)
        for entry in refill:
            self.URL_queue.put(entry)
        return(report)


    def _connect(self,host,user='anonymous',passwd='anonymous@',timeout=20,binary_mode=True):
        ftp=ftplib.FTP(host,user=user,passwd=passwd,timeout=timeout)
        if binary_mode:
            ftp.sendcmd("TYPE i") #swapping out of ASCII mode
        return(ftp)

    def _download_tofile(self,host,file,identity,print_progress=False):
        ftp=self._connect(host)
        lcl_file = (os.path.basename(file))
        #open our write location
        f = open(os.path.join(self.destination,lcl_file),'wb')
        if print_progress:
            print(file)
            self.dl_progress=ProgressBar(ftp.size(file))
            def callback(data):
                f.write(data)
                self.dl_progress.recievebytes(data)
                self.dl_progress.print()
        else:
            def callback(data):
                f.write(data)
        try:
            status = ftp.retrbinary("RETR "+file, callback)
            self.logger.append((identity,lcl_file,status))
        except ftplib.all_errors as e:
            self.logger.append((identity,lcl_file,e))
            status='Failure'
        f.close()
        if '226' not in status:
            os.remove(os.path.join(self.destination,lcl_file))
        ftp.close()



    def Download(self,print_progress=True):
        while not self.URL_queue.empty():
            i,host,file = self.URL_queue.get()
            self._download_tofile(host,file,identity=i,print_progress=print_progress)


def list_dir(address):
    dlist=address.split('//')[-1].split('/')
    host = dlist.pop(0)
    directory = '/'.join(dlist)
    try:
        ftp=ftplib.FTP(host,user='anonymous',passwd='anonymous@',timeout=20)
    except ftplib.all_errors as e:
        raise(e)
    try:
        l=ftp.nlst(directory)
    except ftplib.all_errors as e:
        return(e)
    ftp.close()
    return(l)





































































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
'''