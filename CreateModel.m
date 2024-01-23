
function model=CreateModel()

    % Начальная точка
    xs=6;
    ys=0;
    
    % Целевая точка
    xt=6;
    yt=10;

    % Координаты препятствий
    xobs=[1.5 4.0 1.2 7 5 8 5 7 6];
    yobs=[4.5 3.0 1.5 4 8 6 6 1 2.5];
    robs=[1.5 1.0 0.8 0.8 0.8 0.8 0.8 0.8 0.3];
    
    % Число путевых точек
    n=3;
    
    % Границы пространства
    xmin=00;
    xmax= 10;
    
    ymin=0;
    ymax= 10;

    % Инициализация модели
    model.xs=xs;
    model.ys=ys;
    model.xt=xt;
    model.yt=yt;

    model.xobs=xobs;
    model.yobs=yobs;
    model.robs=robs;

    model.n=n;
    
    model.xmin=xmin;
    model.xmax=xmax;
    model.ymin=ymin;
    model.ymax=ymax;
    
end