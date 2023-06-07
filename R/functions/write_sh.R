write_sh = function(job_name,
                    sh_file,
                    grid_file,
                    inner_file,
                    system = c('sockeye', 'cedar', 'elasti'),
                    env = NA,
                    modules = NA,
                    time = 24, ## in hours
                    mem = 4, ## in GB
                    cpus = 1,
                    gpu = FALSE,
                    partition = FALSE,
                    other_args = NULL
) {
  system = match.arg(system)

  # read the grid
  grid = read.delim(grid_file)

  # set up the script
  if (system == 'sockeye') {
    log_dir = file.path(dirname(base_dir), 'logs', basename(base_dir))
    header_lines = c(
      '#!/bin/bash',
      paste0('#PBS -l walltime=', time, ':00:00,select=1:n',
             ifelse(gpu, 'gpus=', 'cpus='), cpus,
             ':mem=', mem, 'gb'),
      paste0('#PBS -N ', job_name),
      paste0('#PBS -o ', log_dir, '/', job_name, '-^array_index^.out'),
      paste0('#PBS -e ', log_dir, '/', job_name, '-^array_index^.out'),
      ''
    )

    user = system("echo $USER", intern = TRUE)
    env = ifelse(is.na(env), 'CFdb-analysis/env', env)
    env_lines = c(
      '# >>> conda initialize >>>',
      '# !! Contents within this block are managed by \'conda init\' !!',
      paste0('__conda_setup="$(\'/home/', user,
             '/miniconda3/bin/conda\' \'shell.bash\' \'hook\' 2> /dev/null)"'),
      'if [ $? -eq 0 ]; then',
      'eval "$__conda_setup"',
      'else',
      paste0('if [ -f "/home/', user,
             '/miniconda3/etc/profile.d/conda.sh" ]; then'),
      paste0('. "/home/', user, '/miniconda3/etc/profile.d/conda.sh"'),
      'else',
        paste0('export PATH="/home/', user, '/miniconda3/bin:$PATH"'),
      'fi',
      'fi',
      'unset __conda_setup',
      '# <<< conda initialize <<<',
      '',
      paste0('conda activate /project/st-ljfoster-1/', env),
      ''
    )
  } else if (system == 'cedar') {
    log_dir = file.path(dirname(base_dir), 'logs', basename(base_dir))
    header_lines = c(
      '#!/bin/bash',
      paste0('#SBATCH --time=', time, ':00:00'),
      paste0('#SBATCH --job-name=', job_name),
      paste0('#SBATCH --output=', log_dir, '/%x-%j-%a.out'),
      paste0('#SBATCH --mem=', mem, 'G'),
      paste0('#SBATCH --cpus-per-task=', cpus),
      ifelse(partition, paste0('#SBATCH --partition=', partition), ''),
      ifelse(gpu, paste0('#SBATCH --gres=gpu:1'), ''),
      ''
    )

    env_lines = c(
      '',
      if (!is.na(modules))
        paste0(modules, collapse = '\n'),
      '',
      'JOB_SIZE=$1'
    )
  } else if (system == 'elasti') {
    log_dir = file.path('/home/ubuntu/logs/', basename(base_dir))
    header_lines = c(
      '#!/bin/bash',
      paste0('#SBATCH --time=', time, ':00:00'),
      paste0('#SBATCH --job-name=', job_name),
      paste0('#SBATCH --output=', log_dir, '/%x-%j-%a.out'),
      paste0('#SBATCH --mem=', mem, 'G'),
      paste0('#SBATCH --cpus-per-task=', cpus),
      ifelse(gpu, paste0('#SBATCH --gres=gpu:1'), ''),
      ''
    )

    venv = ifelse(is.na(env), '', paste0('source ', env, '/bin/activate'))
    env_lines = switch(gsub("^.*\\.", "", trimws(inner_file)),
           'R' = c(
             '',
             'export PATH=~/R-3.6.0/bin:${PATH}',
             '',
             'JOB_SIZE=$1'
           ),
           'py' = c(
             '',
             venv,
             '',
             'JOB_SIZE=$1'
           )
    )
  } else {
    stop('not sure how to write a sh file for: ', system)
  }

  # set up the final part of the script, which is platform-agnostic
  idx_var = switch(system,
                   'cedar' = 'SLURM_ARRAY_TASK_ID',
                   'elasti' = 'SLURM_ARRAY_TASK_ID',
                   'sockeye' = 'PBS_ARRAY_INDEX')
  run_lines = c(
    'cd ~/git/CFdb-analysis',
    '',
    paste0('START=$((($', idx_var, '-1)*$JOB_SIZE + 1))'),
    paste0('STOP=$((($', idx_var, '-1)*$JOB_SIZE+$JOB_SIZE))'),
    'for i in $( seq $START $STOP ); do',
    paste0('GRID_FILE=', grid_file),
    paste0('LINE_IDX=$((i + 1))'),
    'LINE=`sed "${LINE_IDX}q;d" $GRID_FILE`',
    'IFS=$\'\\t\' PARAMS=($LINE)',
    map_chr(seq_len(ncol(grid)), ~ {
      col_idx = .x
      colname = colnames(grid)[col_idx]
      param_line = paste0(colname, '=${PARAMS[', col_idx - 1, ']}')
      param_line
    }),
    '',
    switch(gsub("^.*\\.", "", trimws(inner_file)),
           'R' = paste0('Rscript ', inner_file, ' \\'),
           'py' = paste0('python ', inner_file, ' \\')
    ),
    map_chr(seq_len(ncol(grid)), ~ {
      col_idx = .x
      colname = colnames(grid)[col_idx]
      if (col_idx < ncol(grid)) {
        arg_line = paste0('  --', colname, ' $', colname, ' \\')
      } else {
        arg_line = paste0('  --', colname, ' $', colname)
      }
      arg_line
    }),
    'done'
  )

  # write to file
  lines = c(header_lines,
            env_lines,
            run_lines)
  sh_file = switch(system,
                   cedar = sh_file,
                   elasti = gsub("\\.sh", "", sh_file) %>%
                       paste0(., '.elasti.sh'),
                   sockeye = gsub("\\.sh", "", sh_file) %>%
                     paste0(., '.torque.sh'))
  sh_dir = dirname(sh_file)
  if (!dir.exists(sh_dir))
    dir.create(sh_dir, recursive = TRUE)
  writeLines(lines, sh_file)
}
