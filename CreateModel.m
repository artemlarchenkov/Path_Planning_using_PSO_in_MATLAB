function model = CreateModel()
    % CREATEMODEL Создает модель для задачи планирования пути
    
    % Начальная точка
    xs = 0;
    ys = 2;
    
    % Целевая точка
    xt = 4;
    yt = 6;
    
    % Координаты препятствий
    xobs = [1.5, 5.0, 1.2];
    yobs = [4.5, 3.0, 1.5];
    
    % Радиусы препятствий
    robs = [1.5, 1.0, 0.8];
    
    % Количество точек управления
    n = 10;
    
    % Границы пространства
    xmin = -10;
    xmax = 10;
    ymin = -10;
    ymax = 10;
    
    % Создание структуры модели
    model.xs = xs;
    model.ys = ys;
    model.xt = xt;
    model.yt = yt;
    model.xobs = xobs;
    model.yobs = yobs;
    model.robs = robs;
    model.n = n;
    model.xmin = xmin;
    model.xmax = xmax;
    model.ymin = ymin;
    model.ymax = ymax;
end
