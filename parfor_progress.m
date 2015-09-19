function percent = parfor_progress(N, caller)
%PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
%   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   PARFOR_PROGRESS(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   PARFOR_PROGRESS updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   PARFOR_PROGRESS(0) deletes parfor_progress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      parfor_progress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress;
%      end
%      parfor_progress(0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/
%Modified by Ruben Pinzon

error(nargchk(0, 2, nargin, 'struct'));

if nargin < 2
    caller = 'Processing';
end

percent = 0;
w = 50; % Width of progress bar

if N > 0
    
    f = fopen('parfor_progress.txt', 'w');
    if f<0
        error('Do you have write permissions for %s?', pwd);
    end
    txt = [caller ' :  0%[>', repmat(' ', 1, w), ']'];
    L = numel(txt); %save the lenght of the text
    if nargout == 0
        disp(txt);
    end
    
    fprintf(f, '%d\t%d\n', N, L); % Save N at the top of progress.txt
    fclose(f);
    
elseif N == 0
    delete('parfor_progress.txt');
    percent = 100;
    
    if nargout == 0
        disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
    end
else
    if ~exist('parfor_progress.txt', 'file')
        error('parfor_progress.txt not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.txt.');
    end
    if nargin < 2
        caller = 'Processing';
    end
    
    f = fopen('parfor_progress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen('parfor_progress.txt', 'r');
    data = fscanf(f, '%d %d', [2 inf]);
    fclose(f);
    
    progress    = data(:,2:end);
    percent     = sum(progress(:))/data(1)*100;
    L           = data(2);
    
    if nargout == 0
        head    = sprintf('%s :%3.0f%%',caller(1:10), percent); % 4 characters wide, percentage
        prog    = sprintf('[%s>%s]',repmat('=', 1, round(percent*w/100)), repmat(' ', 1, w - round(percent*w/100)));
        txt     = [head prog]; 
        
        disp([repmat(char(8), 1, L+2),txt]);
    end
end
