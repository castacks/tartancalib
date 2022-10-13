import numpy as np
import os
import glob

class Experiment:
    def __init__(self,bag,topics,target):
        self.bag = bag 
        self.topics = topics
        self.target = target

data_folder = '/data/final_experiments/'
experiments_output = '/data/final_experiments/'
experiments = []

#ordv single cam
experiments.append(Experiment('ord-6x6-sep12-500f.bag',['/camera_0/image_raw'],'/data/april_6x6.yaml'))
# # experiments.append(Experiment('ord-6x6-sep12-100f.bag',['/camera_0/image_raw'],'/data/april_6x6.yaml'))
experiments.append(Experiment('ord-10x7-sep12-500f.bag',['/camera_0/image_raw'],'/data/400_300.yaml'))


# #ordv multi cam
experiments.append(Experiment('ord-10x7-sep12-500f.bag',['/camera_0/image_raw','/camera_1/image_raw'],'/data/400_300.yaml'))
experiments.append(Experiment('ord-6x6-sep12-500f.bag',['/camera_0/image_raw','/camera_1/image_raw'],'/data/april_6x6.yaml'))

# # # gopro bags
experiments.append(Experiment('gopro-6x6-sep12-500f.bag',['/camera_0/image_raw'],'/data/april_6x6.yaml'))
experiments.append(Experiment('gopro-10x7-sep12-500f.bag',['/camera_0/image_raw'],'/data/400_300.yaml'))

# experiments.append(Experiment('short.bag',['/camera_0/image_raw'],'/data/april_6x6.yaml'))







models = ['ds-none','omni-radtan']

experiment_count = 0
for experiment in experiments:
    for model in models:

        # os.system('rm *.png')
        experiment_folder = experiments_output+'experiment_'+str(experiment_count)
        if (os.path.exists(experiment_folder)):
            filelist = glob.glob(os.path.join(experiment_folder, "*"))
            for f in filelist:
                os.remove(f)
        else:
            os.mkdir(experiment_folder)
        
        command = 'rosrun kalibr tartan_calibrate --bag '+str(data_folder+experiment.bag)+ ' --target ' +  experiment.target+  ' --topics '
        for topic in experiment.topics:
            command += ' '
            command += topic    

        command += ' --models'
        num_cams = len(experiment.topics)
        for _ in range(num_cams):
            command += ' '
            command += model
        
        command += ' --dont-show-report --save_dir '+experiment_folder+'/'
        
        # summary document
        lines = []
        lines.append(experiment.bag)
        lines.append(experiment.target)
        lines.append(str(len(experiment.topics)))
        lines.append(model)
        lines.append('Full command issued: '+command)
        with open(experiment_folder+'/summary.txt','w') as f:
            for line in lines:
                f.write(line)
                f.write('\n')
        print(command)
        os.system(command)


        experiment_count += 1


