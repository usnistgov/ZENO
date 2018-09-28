import os
import glob
import subprocess
import csv

def generate_configs(config_files_dir,
                     bod_files_dir,
                     raw_output_dir,
                     num_tests,
                     num_walks_list):

    bod_file_paths = glob.glob(os.path.join(bod_files_dir, '*.bod'))

    for bod_file_path in sorted(bod_file_paths):

        bod_file_name = os.path.basename(bod_file_path)
        
        bod_prefix = os.path.splitext(bod_file_name)[0]

        seed = 0
        
        for num_walks in num_walks_list:
            
            test_name = bod_prefix + '_' + '%g'%num_walks

            test_input_dir = os.path.join(config_files_dir, test_name)

            test_output_dir = os.path.join(raw_output_dir, test_name)
            
            os.makedirs(test_input_dir, exist_ok=True)

            os.makedirs(test_output_dir, exist_ok=True)
            
            num_interior_samples = num_walks

            for test_num in range(0, num_tests):

                job_name = test_name + '_%03d'%test_num

                config_file_name = job_name + '.cfg'
                
                config_file_path = os.path.join(test_input_dir,
                                                config_file_name)

                csv_file_name = job_name + '.csv'

                csv_file_path = os.path.join(test_output_dir,
                                             csv_file_name)
                
                with open(config_file_path, 'w') as config_file:
                    print('input-file = ' + bod_file_path,
                          file=config_file)
                    print('csv-output-file = ' + csv_file_path,
                          file=config_file)
                    print('num-walks = ' + '%d'%num_walks,
                          file=config_file)
                    print('num-interior-samples = ' + '%d'%num_interior_samples,
                          file=config_file)
                    print('seed = ' + '%d'%seed,
                          file=config_file)
                    print('print-counts',
                          file=config_file)
                    print('print-benchmarks',
                          file=config_file)
        
                seed += 1

def run_tests(config_files_dir,
              zeno_exec):

    for config_file_name in sorted(os.listdir(config_files_dir)):

        config_file_path = os.path.join(config_files_dir,
                                        config_file_name)

        command = [zeno_exec,
                   '-c',
                   config_file_path]

        print(*command, sep=' ')

        subprocess.call(command)

def combine_output(output_csv_files_dir_path,
                   combined_csv_file_path):

    output_csv_file_paths = glob.glob(os.path.join(output_csv_files_dir_path,
                                                   '*.csv'))

    combined_csv_cols = []

    # Add the two header cols
    # Arbitrarily read them from the first output csv file

    combined_csv_cols.append([])
    combined_csv_cols.append([])
        
    with open(output_csv_file_paths[0], newline='') as output_csv_file:
        output_csv_reader = csv.reader(output_csv_file)

        for row in output_csv_reader:
            combined_csv_cols[0].append(row[0])
            combined_csv_cols[1].append(row[1])

    # Add the data col from each output csv file

    for output_csv_file_path in sorted(output_csv_file_paths):
        combined_csv_cols.append([])
        
        with open(output_csv_file_path, newline='') as output_csv_file:
            output_csv_reader = csv.reader(output_csv_file)

            for row in output_csv_reader:
                combined_csv_cols[-1].append(row[2])

    # Write the cols into the combined csv file

    combined_csv_file_dir = os.path.dirname(combined_csv_file_path)

    os.makedirs(combined_csv_file_dir, exist_ok=True)

    combined_csv_rows = [list(i) for i in zip(*combined_csv_cols)]
    
    with open(combined_csv_file_path, 'w', newline='') as combined_csv_file:
        combined_csv_writer = csv.writer(combined_csv_file)

        combined_csv_writer.writerows(combined_csv_rows)
