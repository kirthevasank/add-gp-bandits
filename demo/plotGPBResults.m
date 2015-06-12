% plot GP Bandit results

% Prelims
close all;
PLOT_ERR_BARS = true;
% PLOT_ERR_BARS = false;
NUM_ERR_BARS = 10;
MARKER_SIZE = 8;
LINE_WIDTH = 2;
SAVE_FILE_FORMAT = 'png';

resultsDir = 'results/';
plotColours = {'c', 'b', 'r', 'm', 'k', 'g', [255 128 0]/255, ...
  [76, 0, 153]/253, [102 102 0]/255, 'y'};
plotShapesDot = {'o.', '+.', '*.', 'x.', 's.', 'd.', '^.', 'p.', '>.', 'v.'};
plotShapes = {'o', '+', '*', 'x', 's', 'd', '^', 'p', '>', 'v'};
plotFunc = @semilogx;
plotFunc = @loglog;
plotFunc = @semilogy;
qq = 1:totalNumQueries;
qqq = round(linspace(1,totalNumQueries, NUM_ERR_BARS+2)); qqq = qqq(2:end-1);
numExperiments = sum((randSimpleRegrets(:,1) ~= 0), 1);


% First remove the zero entries
boKDSimpleRegrets = boKDSimpleRegrets(1:numExperiments, :);
boAddSimpleRegrets = boAddSimpleRegrets(1:numExperiments, :, :);
randSimpleRegrets = randSimpleRegrets(1:numExperiments, :);
% Cum Regrets
boKDCumRegrets = boKDCumRegrets(1:numExperiments, :);
boAddCumRegrets = boAddCumRegrets(1:numExperiments, :, :);
randCumRegrets = randCumRegrets(1:numExperiments, :);

for regIter = 1:2

  if regIter == 1 % First Do Simple Regret
    figTitlePrefix = 'Simple-Regret';
    % Mean
    KDRegMean = mean(boKDSimpleRegrets, 1);
    AddRegMean = mean(boAddSimpleRegrets, 1);
    eiRegMean = mean(boEISimpleRegrets, 1);
    randRegMean = mean(randSimpleRegrets, 1);
    % Std
    KDRegStdErr = std(boKDSimpleRegrets, 1)/sqrt(numExperiments);
    AddRegStdErr = std(boAddSimpleRegrets, 1)/sqrt(numExperiments);
    eiRegStdErr = std(boEISimpleRegrets, 1)/sqrt(numExperiments);
    randRegStdErr = std(randSimpleRegrets, 1)/sqrt(numExperiments);
    % For diRect
    diRectReg = diRectSimpleRegret;

  else % Now do Cumulative Regret
    figTitlePrefix = 'Cumulative-Regret';
    % Mean
    KDRegMean = mean(boKDCumRegrets, 1);
    AddRegMean = mean(boAddCumRegrets, 1);
    eiRegMean = mean(boEICumRegrets, 1);
    randRegMean = mean(randCumRegrets, 1);
    % Std
    KDRegStdErr = std(boKDCumRegrets, 1)/sqrt(numExperiments);
    AddRegStdErr = std(boAddCumRegrets, 1)/sqrt(numExperiments);
    eiRegStdErr = std(boEIAddRegrets, 1)/sqrt(numExperiments);
    randRegStdErr = std(randCumRegrets, 1)/sqrt(numExperiments);
    % For diRect
    diRectReg = diRectCumRegret;

  end

  % Correct extra terms for diRect
  diRectReg = diRectReg(1:totalNumQueries);

  % Plot Iteration statistics
  figure;
  plotFunc(qqq, randRegMean(qqq), plotShapes{3}, 'Color', plotColours{3}, ...
     'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH); hold on,
  if regIter == 1, % Don't plot Cum Regret for DiRect
    plotFunc(qqq, diRectReg(qqq), plotShapes{4}, 'Color', plotColours{4}, ...
      'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH); hold on,
  end
  plotFunc(qqq, eiRegMean(qqq), plotShapes{2}, 'Color', plotColours{2}, ...
    'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH); hold on,
  plotFunc(qqq, KDRegMean(qqq), plotShapes{1}, 'Color', plotColours{1}, ...
    'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH); hold on,
  if regIter == 1
    legEntries = {'Random', 'DiRect', 'BO-EI', 'BO-KD'};
  else
    legEntries = {'Random', 'BO-EI', 'BO-KD'};
  end
  numBaseLegEntries = numel(legEntries);
  for i = 1:numdCands
    plotFunc(qqq, AddRegMean(1,qqq,i), plotShapes{4+i}, 'Color', plotColours{4+i}, ...
  'MarkerSize', MARKER_SIZE, 'LineWidth', LINE_WIDTH);
    legEntries{numBaseLegEntries+i} = sprintf('BO-%d', numDimsPerGroupCands(i));
  end
  legend(legEntries);

  % Now reproduce the curve without the bullets
  plotFunc(qq, randRegMean, 'Color', plotColours{3}, 'LineWidth', LINE_WIDTH); hold on,
  if regIter == 1, % Don't plot Cum Regret for DiRect
    plotFunc(qq, diRectReg, 'Color', plotColours{4}, 'LineWidth', LINE_WIDTH); hold on,
  end
  plotFunc(qq, eiRegMean, 'Color', plotColours{2}, 'LineWidth', LINE_WIDTH); hold on,
  plotFunc(qq, KDRegMean, 'Color', plotColours{1}, 'LineWidth', LINE_WIDTH); hold on,
  for i = 1:numdCands
    plotFunc(qq, AddRegMean(1,:,i), 'Color', plotColours{4+i}, 'LineWidth', LINE_WIDTH);
  end

  % Plot Error Bars
  if PLOT_ERR_BARS & (numExperiments > 1)
    errorbar(qqq, randRegMean(qqq), randRegStdErr(qqq), '.', 'Color', plotColours{3});
    errorbar(qqq, eiRegMean(qqq), eiRegStdErr(qqq), '.', 'Color', plotColours{2});
    errorbar(qqq, KDRegMean(qqq), KDRegStdErr(qqq), '.', 'Color', plotColours{1});
    for i = 1:numdCands
      errorbar(qqq, AddRegMean(1,qqq,i), AddRegStdErr(1,qqq,i), '.', ...
        'Color', plotColours{4+i});
    end
  end

  % the minimum and maximum for Plotting
  addRegMeanMaxVals = AddRegMean(1, 1, :);
  maxPlotVal = max([KDRegMean(1); randRegMean(1); eiRegMean(1); ...
    addRegMeanMaxVals(:)]);
  addRegMeanMinVals = AddRegMean(1, end, :);
  minPlotVal = min([KDRegMean(end); randRegMean(end); eiRegMean(end); ...
    addRegMeanMinVals(:)]);
  
  plotRange = maxPlotVal - minPlotVal; 
  xlim([0 1.05*totalNumQueries]);
  ylim([minPlotVal - plotRange*0.02, maxPlotVal + plotRange*0.5]);
  ylabel(figTitlePrefix);
  xlabel('Number of Queries (T)');
  titlestr = sprintf('(D,d'',M'') = (%d,%d,%d)', numDims, ...
    trueNumDimsPerGroup, floor(numDims/trueNumDimsPerGroup));
  title(titlestr);

  titleStr = sprintf('%s, D = %d', figTitlePrefix, numDims);
  title(titleStr);


end
