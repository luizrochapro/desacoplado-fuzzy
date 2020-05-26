#coding: utf-8
class Log:
    
    log_file = None
    path_file = None

    
    def __init__(self, path_file):
        self.path_file  = path_file

    def open_file(self):
        self.log_file = open(self.path_file,'w', encoding='utf-8')
    
    def close_file(self):
        self.log_file.close()

    def write_log(self, str):
        self.log_file.write(str)
        self.log_file.write('\n')
    
    def write_log_space(self):
        self.log_file.write('\n')
    
    def write_log_div(self):
        self.log_file.write('=============================================================================================================' + '\n')

    def write_log_iter_mark(self, iter):
        self.log_file.write('=============================================================================================================' + '\n')
        self.log_file.write('==========================       ITERAÇÃO  {:4}        ======================================================'.format(str(iter)) + '\n')
        self.log_file.write('=============================================================================================================' + '\n')

    def write_log_list_fuzzy(self, lista):
        for item in lista:
            self.log_file.write(str(item.f)) 
            self.log_file.write('\n')       