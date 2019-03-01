

[pValues, updateValues, C] = fastpolarize(10,1);
map = polarizeExample(10,1,...
    'pValues',pValues,'updateValues',updateValues);
saveas(gcf,'polargif/0.jpg')

for i =1:30
    map = polarizeExample(10,min(max(i,0),1),'initialMap',map,...
        'pValues',pValues,'updateValues',updateValues);
    saveas(gcf,strcat('polargif/', num2str(i), '.jpg'))
end

close all
