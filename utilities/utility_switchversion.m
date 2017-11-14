% a q&d routine to switch between tslbatch routines
rootp='d:';

% the ones to kill
rmpath([rootp '\hh\matlab\beastybites\tslbatch']);
rmpath([rootp '\hh\matlab\beastybites\tslbatch\func']);
rmpath([rootp '\hh\matlab\beastybites\tslbatch\utilities']);

% these ones instead
addpath([rootp '\hh\matlab\beastybites\tslbatch\_archive\V1.1\']);
addpath([rootp '\hh\matlab\beastybites\tslbatch\_archive\V1.1\func']);
addpath([rootp '\hh\matlab\beastybites\tslbatch\_archive\V1.1\utilities']);

