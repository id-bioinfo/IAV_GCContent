clear all

colors(3,:) = [105,105,105]/255; % dim grey

filelist = {
'svm_Hu1_year.txt';
'svm_Hu2_year.txt';
'svm_Hu3_year.txt';
'svm_Hu4_year.txt';
'svm_Ca1_year.txt';
'svm_Ca2_year.txt';
'svm_Eq_year.txt';
'svm_Sw1_year.txt';
'svm_Sw2_year.txt';
'svm_Sw3_year.txt';
'svm_Sw4_year.txt';
'svm_Sw5_year.txt';
 };

prefixPath = '\svm_result\';

figure(1)
% ymin = 0.38;
% ymax = 0.45;
slope = zeros(size(filelist,1),1);
pvalue = zeros(size(filelist,1),1);
bints = zeros(size(filelist,1),2);
R2 =  zeros(size(filelist,1),1);

for i=1:size(filelist,1)
    filepath = strcat(prefixPath,filelist{i,1});

    d = load(filepath);
  
     if i == 1
         d = d(find(d(:,1)<2010), :);
     end
    mysize = int2str(size(d,1));
        

    x = d(:,1);
    y = d(:,2);
    
    subplot(2, 6, i);

    plot(x,y,'.', 'Color',colors(3,:));
    hold on

     titlename = strrep(filelist{i,1},'_','/');
     titlename = strrep(titlename,'/year','');
     titlename = strrep(titlename,'svm/','');
     titlename = strrep(titlename,'.txt','');
 
     titlename = strcat(titlename, ' (', mysize, ')');
     title(titlename, 'FontSize', 7);  
     
      xmin=floor(min(x(:,1)));
      xmax=ceil(max(x(:,1)));
        
      xlim([xmin,xmax]);
      interval =  (xmax-xmin)/2;
      set(gca,'XTick',[xmin:interval:xmax]);
    set(gca,'ylim', [-1, 1], 'FontSize', 7);  
    
     if i==1
         ylabel('Classification', 'FontSize', 7);  
     end
end

%print('Classification_predict.pdf','-dpdf', '-r300');

