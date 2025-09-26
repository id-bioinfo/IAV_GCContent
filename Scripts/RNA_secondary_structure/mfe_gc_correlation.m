clear all

colors(3,:) = [105,105,105]/255; % dim grey

filelist = {
'_rnafold_addGC3_12_statics_Hu1.csv';
'_rnafold_addGC3_12_statics_Hu2.csv';
'_rnafold_addGC3_12_statics_Hu3.csv';
'_rnafold_addGC3_12_statics_Hu4.csv';
'_rnafold_addGC3_12_statics_Eq.csv';
'_rnafold_addGC3_12_statics_Ca1.csv';
'_rnafold_addGC3_12_statics_Ca2.csv';
'_rnafold_addGC3_12_statics_Sw1.csv';
'_rnafold_addGC3_12_statics_Sw2.csv';
'_rnafold_addGC3_12_statics_Sw3.csv';
'_rnafold_addGC3_12_statics_Sw4.csv';
'_rnafold_addGC3_12_statics_Sw5.csv';
%'_rnafold_statics_Sp.csv';
%'_rnafold_statics_Av.csv';
 };

prefixPath = 'MFE/NS1/NS1';

 % NS
 ymin = -280;
 ymax = -200;
 xmin = 40;
 xmax = 45;

slope = zeros(size(filelist,1),1);
pvalue = zeros(size(filelist,1),1);
bints = zeros(size(filelist,1),2);
R2 =  zeros(size(filelist,1),1);
corr_value = zeros(size(filelist,1),1);
corr_pvalue = zeros(size(filelist,1),1);
counts = zeros(size(filelist,1),1);
gc_content = zeros(size(filelist,1),1);
gc3_content = zeros(size(filelist,1),1);
gc3_corr_value = zeros(size(filelist,1),1);
gc3_corr_pvalue = zeros(size(filelist,1),1);
gc12_content = zeros(size(filelist,1),1);
gc12_corr_value = zeros(size(filelist,1),1);
gc12_corr_pvalue = zeros(size(filelist,1),1);
mfe = zeros(size(filelist,1),1);
std_gc = zeros(size(filelist,1),1);
std_mfe = zeros(size(filelist,1),1);

for i=1:size(filelist,1)
    filepath = strcat(prefixPath,filelist{i,1});
    alld = importdata(filepath);
    d = alld.data;
    %d = load(filepath);
  
     if i == 1
         d = d(find(d(:,1)<2010), :);
     end
       
     %d = d(find(d(:,5)==890), :);
     
     if(size(d,1) < 1) 
         continue;
     end
    mysize = int2str(size(d,1));
    
    x = d(:,3); %GC
    %x = d(:,1); %time
    y = d(:,2); %GC
    %y = d(:,2); %MFE

    x_gc3 = d(:,6); %GC
    x_gc12 = d(:,7); %GC

    avgx = x;
    avgy = y;

    p=polyfit(avgx,avgy,1);  
    X = [ones(size(avgx)), avgx];
    [b,bint,r,rint,stats] = regress(avgy, X);
        
    x1=linspace(min(x),max(x));  
    y1=polyval(p,x1);  
    slope(i,1) = p(1,1);
    pvalue(i,1) = stats(1,3);
    bints(i,1) = bint(2,1);
    bints(i,2) = bint(2,2);
    R2(i,1) = stats(1,1);
    [rho,pval] = corr(avgx,avgy);
    corr_value(i,1) = rho;
    corr_pvalue(i,1) = pval;
    counts(i,1) = size(d,1);
    gc_content(i,1) = mean(x);
    mfe(i,1) = mean(y);
    std_gc(i,1) = std(x);
    std_mfe(i,1) = std(y);

    [rho,pval] = corr(x_gc3,avgy);
    gc3_corr_value(i,1) = rho;
    gc3_corr_pvalue(i,1) = pval;

    [rho,pval] = corr(x_gc12,avgy);
    gc12_corr_value(i,1) = rho;
    gc12_corr_pvalue(i,1) = pval;


    %subplot(3, 4, i);

    %plot(x,y,'.', 'Color',colors(3,:));
    hold on
     if i <= 12 && i >= 1
     %plot(x1,y1, 'r');
     end
     
     if i == 15
        d1 = d(find(d(:,1)<=1958), :);
        %d1 = d1(find(d1(:,1)>1918), :);
        avgx = d1(:,3);
        avgy = d1(:,2);
        p=polyfit(avgx,avgy,1);  
        X = [ones(size(avgx)), avgx];
        [b,bint,r,rint,stats] = regress(avgy, X);
        
        x1=linspace(min( avgx),max(avgx));  
        y1=polyval(p,x1);  
        slope(i,1) = p(1,1);
        pvalue(i,1) = stats(1,3);
		bints(i,1) = bint(2,1);
		bints(i,2) = bint(2,2);
		R2(i,1) = stats(1,1);
      
		plot(x1,y1, 'r');
      
		hold on
      
		d2 = d(find(d(:,1)>1958), :);
        avgx = d2(:,3);
        avgy = d2(:,2);%*100;
        p=polyfit(avgx,avgy,1);  
        X = [ones(size(avgx)), avgx];
        [b,bint,r,rint,stats] = regress(avgy, X);
        
        x1=linspace(min( avgx),max(avgx));  
		y1=polyval(p,x1);  
		slope(i,2) = p(1,1);
		pvalue(i,2) = stats(1,3);
		bints(i,3) = bint(2,1);
		bints(i,4) = bint(2,2);
		R2(i,2) = stats(1,1);
      
		%plot(x1,y1, 'r');

    end
     
    titlename = strrep(filelist{i,1},'_rnafold_statics_','');
    titlename = strrep(titlename,'.csv','');
 
    titlename = strcat(titlename, ' (', mysize, ')');
    title(titlename, 'FontSize', 7);  
    
%     xmin=floor(min(x(:,1)));
%       xmax=ceil(max(x(:,1)));
%         
%       xlim([xmin,xmax]);
%       interval =  (xmax-xmin)/2;
%       set(gca,'XTick',[xmin:interval:xmax]);

    set(gca,'xlim', [xmin, xmax], 'FontSize', 6);
    set(gca,'ylim', [ymin, ymax], 'FontSize', 6);
     
     if i==1
        %ylabel('GC content (%)', 'FontSize', 7);  
        ylabel('MFE (kcal/mol)', 'FontSize', 7);  
		xlabel('GC content (%)', 'FontSize', 7);  
        %xlabel('Time', 'FontSize', 7);  
     end
end

summary = zeros(size(filelist,1),5);
summary(:,1) = slope(:,1);
summary(:,2:3) = bints(:,1:2);
summary(:,4) = pvalue(:,1);
summary(:,5) =R2(:,1);
summary(:,6) = corr_value(:,1);
summary(:,7) = corr_pvalue(:,1);
summary(:,8) = counts(:,1);
summary(:,9) = gc_content(:,1);
summary(:,10) = mfe(:,1);
summary(:,11) = std_gc(:,1);
summary(:,12) = std_mfe(:,1);
summary(:,13) = gc3_corr_value(:,1);
summary(:,14) = gc3_corr_pvalue(:,1);
summary(:,15) = gc12_corr_value(:,1);
summary(:,16) = gc12_corr_pvalue(:,1);


