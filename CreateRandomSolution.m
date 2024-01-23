function sol1 = CreateRandomSolution(model)
    % CREATERANDOMSOLUTION Создает случайное решение для задачи планирования пути
    
    % Извлечение параметров модели
    n = model.n;
    xmin = model.xmin;
    xmax = model.xmax;
    ymin = model.ymin;
    ymax = model.ymax;

    % Генерация случайных координат для точек пути
    sol1.x = unifrnd(xmin, xmax, 1, n);
    sol1.y = unifrnd(ymin, ymax, 1, n);
end
