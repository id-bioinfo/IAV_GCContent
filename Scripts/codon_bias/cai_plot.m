clear all

colors(3,:) = [105,105,105]/255; % dim grey

filelist = {
'Hu1_cai_time.csv';
'Hu2_cai_time.csv';
'Hu3_cai_time.csv';
'Hu4_cai_time.csv';
'Eq_cai_time.csv';
'Ca1_cai_time.csv';
'Ca2_cai_time.csv';
'Sw1_cai_time.csv';
'Sw2_cai_time.csv';
'Sw3_cai_time.csv';
'Sw4_cai_time.csv';
'Sw5_cai_time.csv';
 };

prefixPath = '../../CAI/';

figure(1)
% ymin = 0.38;
% ymax = 0.45;
slope = zeros(size(filelist,1),1);
pvalue = zeros(size(filelist,1),1);
bints = zeros(size(filelist,1),2);
R2 =  zeros(size(filelist,1),1);

for i=1:size(filelist,1)
    filepath = strcat(prefixPath,filelist{i,1});
    alld = importdata(filepath);
    d = alld.data;
    %d = load(filepath);
  
     if i == 1
         d = d(find(d(:,2)<2010), :);
     end
    mysize = int2str(size(d,1));


    x = d(:,2);
    y = d(:,1);

    p=polyfit(x,y,1);  
    X = [ones(size(x)), x];
    [b,bint,r,rint,stats] = regress(y, X);
        
    x1=linspace(min(x),max(x));  
    y1=polyval(p,x1);  
    slope(i,1) = p(1,1);
    pvalue(i,1) = stats(1,3);
    bints(i,1) = bint(2,1);
    bints(i,2) = bint(2,2);
    R2(i,1) = stats(1,1);
    
    subplot(3, 4, i);

    plot(x,y,'.', 'Color',colors(3,:));
    hold on

     titlename = strrep(filelist{i,1},'_','');
     titlename = strrep(titlename,'time','');
     titlename = strrep(titlename,'cai','');
     titlename = strrep(titlename,'.csv','');
 
     titlename = strcat(titlename, ' (', mysize, ')');
     title(titlename, 'FontSize', 7);  
     
      xmin=floor(min(x(:,1)));
      xmax=ceil(max(x(:,1)));
        
      xlim([xmin,xmax]);
      if i <=4
        ylim([0.7,0.73]);
      elseif i <= 7
         ylim([0.62,0.67]);
      else
         ylim([0.53,0.57]); 
      end

      interval =  (xmax-xmin)/2;
      set(gca,'XTick',[xmin:interval:xmax]);

     if i==1
         ylabel('CAI', 'FontSize', 7);  
         xlabel('time', 'FontSize', 6);  
     end
end

print('cai_lineage.pdf','-dpdf', '-r300');

