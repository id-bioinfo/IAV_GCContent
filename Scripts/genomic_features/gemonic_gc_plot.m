clear all

colors(3,:) = [105,105,105]/255; % dim grey

filelist = {
'Hu1';
'Hu2';
'Hu3';
'Hu4';
'Eq';
'Ca1';
'Ca2';
'Sw1';
'Sw2';
'Sw3';
'Sw4';
'Sw5';
 };

prefixPath = '../../GC_content/';
figure(1)

slope = zeros(size(filelist,1),1);
pvalue = zeros(size(filelist,1),1);
bints = zeros(size(filelist,1),2);
R2 =  zeros(size(filelist,1),1);

for i=1:size(filelist,1)
    filepath = strcat(prefixPath,filelist{i,1},'_gc_genomic_year_name.txt');
    %rd = importdata(filepath);
    rd = readtable(filepath, 'FileType', 'text', 'Delimiter', '\t');
    d = table2array(rd(:,1:2));
    %d = load(filepath)
  
    if i == 1
        d = d(find(d(:,1)<2010), :);
    end
    mysize = int2str(size(d,1));

    x = d(:,1);
    y = d(:,2);
    
    p=polyfit(x,y,1);  
    X = [ones(size(x)), y];
    [b,bint,r,rint,stats] = regress(y, X);
        
    x1=linspace(min(x),max(x));  
    y1=polyval(p,x1);  
    slope(i,1) = p(1,1);
    pvalue(i,1) = stats(1,3);
    bints(i,1) = bint(2,1);
    bints(i,2) = bint(2,2);
    R2(i,1) = stats(1,1);
    
    subplot(3,4, i);

    plot(x,y,'.', 'Color',colors(3,:));

    hold on
    if i <= 12 && i > 1
        plot(x1,y1, 'r');
    end
     
    if i == 1
        d1 = d(find(d(:,1)<=1957), :);
        x1 = d1(:,1);
        y1 = d1(:,2);
        p=polyfit(x1,y1,1);  
        X = [ones(size(x1)), x1];
        [b,bint,r,rint,stats] = regress(y1, X);
        
        x1=linspace(min(x1),max(x1));  
        y1=polyval(p,x1);  
        slope(i,1) = p(1,1);
        pvalue(i,1) = stats(1,3);
        bints(i,1) = bint(2,1);
        bints(i,2) = bint(2,2);
        R2(i,1) = stats(1,1);
      
        plot(x1,y1, 'r');      
        hold on
        
        d2 = d(find(d(:,1)>=1977), :);
        x2 = d2(:,1);
        y2 = d2(:,2);
        p=polyfit(x2,y2,1);  
        X = [ones(size(x2)), x2];
        [b,bint,r,rint,stats] = regress(y2, X);
        
        x2=linspace(min( x2),max(x2));  
        y2=polyval(p,x2);  
        slope(i,2) = p(1,1);
        pvalue(i,2) = stats(1,3);
        bints(i,3) = bint(2,1);
        bints(i,4) = bint(2,2);
        R2(i,2) = stats(1,1);
      
        plot(x2,y2, 'r');
    
     end
     
    titlename = strrep(filelist{i,1},'_year.txt',''); 
    titlename = strcat(titlename, ' (', mysize, ')');
    title(titlename, 'FontSize', 7);  
     
    xmin=floor(min(x(:,1)));
    xmax=ceil(max(x(:,1)));
        
    xlim([xmin,xmax]);
    interval =  (xmax-xmin)/2;
    set(gca,'XTick',[xmin:interval:xmax]);

    set(gca,'ylim', [41.5, 46], 'FontSize', 7);  %genomic GC
     
    if i==1
        ylabel('Genomic GC Content (%)', 'FontSize', 7);  
    end
end

summary = zeros(13,5);
summary(1,1) = slope(1,1);
summary(1,2:3) = bints(1,1:2);
summary(1,4) = pvalue(1,1);
summary(1,5) =R2(1,1);

summary(2,1) = slope(1,2);
summary(2,2:3) = bints(1,3:4);
summary(2,4) = pvalue(1,2);
summary(2,5) =R2(1,2);

summary(3:13,1) = slope(2:12,1);
summary(3:13,2:3) = bints(2:12,1:2);
summary(3:13,4) = pvalue(2:12,1);
summary(3:13,5) =R2(2:12,1);

print('Genomic_GC_Content_2Hu1.pdf','-dpdf', '-r300');
