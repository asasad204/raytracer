// ======================================================================
//  Core Classes and Constants
// ======================================================================
const EPSILON = 0.001;

class Vec {
    constructor(x, y, z) {
        this.x = x || 0;
        this.y = y || 0;
        this.z = z || 0;
    }

    add(v) { 
        if (!v) return this;
        return new Vec(this.x + (v.x || 0), this.y + (v.y || 0), this.z + (v.z || 0)); 
    }
    
    sub(v) { 
        if (!v) return this;
        return new Vec(this.x - (v.x || 0), this.y - (v.y || 0), this.z - (v.z || 0)); 
    }
    
    mul(s) {
        if (s instanceof Vec) {
            return new Vec(this.x * (s.x || 0), this.y * (s.y || 0), this.z * (s.z || 0));
        }
        return new Vec(this.x * (s || 0), this.y * (s || 0), this.z * (s || 0));
    }
    
    dot(v) { 
        if (!v) return 0;
        return this.x * (v.x || 0) + this.y * (v.y || 0) + this.z * (v.z || 0); 
    }
    
    cross(v) {
        if (!v) return new Vec();
        return new Vec(
            this.y * (v.z || 0) - this.z * (v.y || 0),
            this.z * (v.x || 0) - this.x * (v.z || 0),
            this.x * (v.y || 0) - this.y * (v.x || 0)
        );
    }
    
    length() { 
        return Math.sqrt(this.dot(this)); 
    }
    
    normalized() { 
        const len = this.length();
        if (len <= 0) return new Vec();
        return this.mul(1 / len); 
    }
}

function ReflectRay(v, n) {
    if (!v || !n) return new Vec();
    return n.mul(2 * v.dot(n)).sub(v);
}

class Color {
    constructor(r, g, b) {
        this.r = Math.max(0, Math.min(255, Number(r) || 0));
        this.g = Math.max(0, Math.min(255, Number(g) || 0));
        this.b = Math.max(0, Math.min(255, Number(b) || 0));
    }

    mul(s) {
        if (!s && s !== 0) return new Color(0, 0, 0);
        
        if (s instanceof Color) {
            return new Color(
                (this.r * s.r) / 255,
                (this.g * s.g) / 255,
                (this.b * s.b) / 255
            );
        }
        
        const scalar = Number(s);
        if (isNaN(scalar)) return new Color(0, 0, 0);
        
        return new Color(
            Math.min(255, this.r * scalar),
            Math.min(255, this.g * scalar),
            Math.min(255, this.b * scalar)
        );
    }

    add(c) {
        if (!c) return new Color(this.r, this.g, this.b);
        return new Color(
            Math.min(255, this.r + (c.r || 0)),
            Math.min(255, this.g + (c.g || 0)),
            Math.min(255, this.b + (c.b || 0))
        );
    }
}

// ======================================================================
//  Scene Objects
// ======================================================================

class Triangle {
    constructor(v0, v1, v2, color = null) {
        this.vertices = [v0, v1, v2];
        this.type = 'triangle';
        this.normal = this.calculateNormal();
        
        // Assign a random color if none provided
        this.color = color || new Color(
            Math.floor(Math.random() * 256),
            Math.floor(Math.random() * 256),
            Math.floor(Math.random() * 256)
        );
        
        this.material = {
            diffuse: this.color,
            specular: 500,
            reflective: 0.2
        };
    }

    calculateNormal() {
        const edge1 = this.vertices[1].sub(this.vertices[0]);
        const edge2 = this.vertices[2].sub(this.vertices[0]);
        return edge1.cross(edge2).normalized();
    }
}

class Sphere {
    constructor(center, radius, color, specular, reflective) {
        this.type = 'sphere';
        this.center = center;
        this.radius = radius;
        this.color = color;
        this.specular = specular;
        this.reflective = reflective;
        this.material = {
            diffuse: color,
            specular: specular,
            reflective: reflective
        };
    }
}

// ======================================================================
//  Scene Configuration
// ======================================================================

const LightType = {
    AMBIENT: 0,
    POINT: 1,
    DIRECTIONAL: 2
};

const scene = {
    viewport_size: 1,
    projection_plane_z: 1,
    camera_position: new Vec(0, 0, 0),
    background_color: new Color(0, 0, 0),
    recursion_limit: 3,
    meshes: [],
    lights: [
        { type: LightType.AMBIENT, intensity: 0.2 },
        { type: LightType.POINT, intensity: 0.6, position: new Vec(2, 1, 0) },
        { type: LightType.DIRECTIONAL, intensity: 0.2, position: new Vec(1, 4, 4) }
    ],
    spheres: [
        new Sphere(new Vec(0, -1, 3), 1, new Color(255, 0, 0), 500, 0.2),
        new Sphere(new Vec(-2, 0, 4), 1, new Color(0, 255, 0), 10, 0.4),
        new Sphere(new Vec(2, 0, 4), 1, new Color(0, 0, 255), 500, 0.3),
        new Sphere(new Vec(0, -5001, 0), 5000, new Color(255, 255, 0), 1000, 0.5)
    ]
};

// ======================================================================
//  Canvas Rendering
// ======================================================================

const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');
let canvas_buffer = ctx.getImageData(0, 0, canvas.width, canvas.height);

function putPixel(x, y, color) {
    x = canvas.width / 2 + Math.floor(x);
    y = canvas.height / 2 - Math.floor(y) - 1;

    if (x < 0 || x >= canvas.width || y < 0 || y >= canvas.height) return;

    const offset = 4 * (x + canvas_buffer.width * y);
    canvas_buffer.data[offset] = color.r;
    canvas_buffer.data[offset + 1] = color.g;
    canvas_buffer.data[offset + 2] = color.b;
    canvas_buffer.data[offset + 3] = 255;
}

function updateCanvas() {
    ctx.putImageData(canvas_buffer, 0, 0);
}

function clearCanvas() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    canvas_buffer = ctx.createImageData(canvas.width, canvas.height);
}

// ======================================================================
//  Ray Intersection Functions
// ======================================================================
function intersectRaySphere(origin, direction, sphere) {
    const CO = origin.sub(sphere.center);
    const a = direction.dot(direction);
    const b = 2 * CO.dot(direction);
    const c = CO.dot(CO) - sphere.radius * sphere.radius;
    const discriminant = b * b - 4 * a * c;

    if (discriminant < 0) return null;

    const sqrtDisc = Math.sqrt(discriminant);
    const t1 = (-b + sqrtDisc) / (2 * a);
    const t2 = (-b - sqrtDisc) / (2 * a);

    if (t2 > EPSILON) return t2;
    if (t1 > EPSILON) return t1;

    return null;
}

function intersectRayTriangle(origin, direction, triangle) {
    const [v0, v1, v2] = triangle.vertices;
    const edge1 = v1.sub(v0);
    const edge2 = v2.sub(v0);
    const h = direction.cross(edge2);
    const a = edge1.dot(h);

    if (a > -EPSILON && a < EPSILON) return null;

    const f = 1 / a;
    const s = origin.sub(v0);
    const u = f * s.dot(h);

    if (u < 0 || u > 1) return null;

    const q = s.cross(edge1);
    const v = f * direction.dot(q);

    if (v < 0 || u + v > 1) return null;

    const t = f * edge2.dot(q);
    if (t > EPSILON) return t;

    return null;
}

// ======================================================================
//  OBJ Loading
// ======================================================================

function parseOBJ(text) {
    const vertices = [];
    const faces = [];

    const lines = text.split('\n');
    for (const line of lines) {
        const trimmed = line.trim();
        if (!trimmed || trimmed.startsWith('#')) continue;

        const tokens = trimmed.split(/\s+/);
        const keyword = tokens[0];
        const args = tokens.slice(1);

        switch (keyword) {
            case 'v': // Vertex
                vertices.push(new Vec(
                    parseFloat(args[0]),
                    parseFloat(args[1]),
                    parseFloat(args[2])
                ));
                break;
            case 'f': // Face
                const faceVertices = [];
                for (const arg of args) {
                    const vertexIndex = parseInt(arg.split('/')[0]) - 1;
                    if (!isNaN(vertexIndex)) {
                        faceVertices.push(vertexIndex);
                    }
                }
                faces.push(faceVertices);
                break;
        }
    }

    return { vertices, faces };
}

async function loadScene() {
    const objFile = document.getElementById('objFile').files[0];
    if (!objFile) {
        alert('Please select an OBJ file');
        return;
    }

    try {
        const objText = await objFile.text();
        const objData = parseOBJ(objText);

        scene.meshes = [];
        scene.spheres = [];
        for (const face of objData.faces) {
            for (let i = 2; i < face.length; i++) {
                const v0 = objData.vertices[face[0]];
                const v1 = objData.vertices[face[i-1]];
                const v2 = objData.vertices[face[i]];
                
                if (v0 && v1 && v2) {
                    scene.meshes.push(new Triangle(v0, v1, v2, new Color(255, 255, 255)));
                }
            }
        }
        
        // Auto-adjust camera
        autoPositionCamera(objData.vertices);
        render();
    } catch (error) {
        console.error('Error loading scene:', error);
        alert('Error loading scene: ' + error.message);
    }
}

function autoPositionCamera(vertices) {
    if (vertices.length === 0) return;
    
    let min = new Vec(Infinity, Infinity, Infinity);
    let max = new Vec(-Infinity, -Infinity, -Infinity);
    
    vertices.forEach(v => {
        min = new Vec(Math.min(min.x, v.x), Math.min(min.y, v.y), Math.min(min.z, v.z));
        max = new Vec(Math.max(max.x, v.x), Math.max(max.y, v.y), Math.max(max.z, v.z));
    });
    
    const center = new Vec(
        (min.x + max.x) / 2,
        (min.y + max.y) / 2,
        (min.z + max.z) / 2
    );
    
    const size = Math.max(max.x - min.x, max.y - min.y, max.z - min.z);
    
    scene.camera_position = new Vec(
        center.x,
        center.y,
        center.z - size * 2
    );
    
    scene.viewport_size = size/4; // Adjust viewport size based on object size

}

// ======================================================================
//  Main Rendering Function with Web Workers
// ======================================================================

function render() {
    clearCanvas();
    scene.recursion_limit = parseInt(document.getElementById('rec-limit').value) || 3;

    const width = canvas.width;
    const height = canvas.height;
    const numWorkers = navigator.hardwareConcurrency || 4;
    console.log(`Using ${numWorkers} workers for rendering`);
    const rowsPerWorker = Math.ceil(height / numWorkers);
    let completedWorkers = 0;

    // Create worker code as a string
    const workerCode = `
        const EPSILON = ${EPSILON};

        ${Vec.toString()}
        ${Color.toString()}
        ${ReflectRay.toString()}
        ${intersectRayTriangle.toString()}
        ${intersectRaySphere.toString()}

        function traceRay(origin, direction, min_t, max_t, depth, scene) {
            const intersection = closestIntersection(origin, direction, min_t, max_t, scene);
            if (!intersection) return scene.background_color;

            const { obj, t } = intersection;
            const point = origin.add(direction.mul(t));

            let normal;
            if (obj.type === 'triangle') {
                const faceNormal = obj.normal.normalized();
                normal = (faceNormal.dot(direction)) > 0 ? faceNormal.mul(-1) : faceNormal;
            } else if (obj.type === 'sphere') {
                normal = point.sub(obj.center).normalized();
            }

            const view = direction.mul(-1);
            const lighting = computeLighting(point, normal, view, obj.material.specular, scene);
            const local_color = obj.material.diffuse.mul(lighting);

            if (depth <= 0 || !obj.material.reflective) {
                return local_color;
            }

            const reflected_ray = ReflectRay(view, normal);
            const reflected_color = traceRay(point, reflected_ray, EPSILON, Infinity, depth - 1, scene);
            
            return local_color.mul(1 - obj.material.reflective)
                        .add(reflected_color.mul(obj.material.reflective));
        }

        function canvasToViewport(x, y, width, height, scene) {
            return new Vec(
                x * scene.viewport_size / width,
                y * scene.viewport_size / height,
                scene.projection_plane_z
            );
        }

        function closestIntersection(origin, direction, min_t, max_t, scene) {
            let closest_t = Infinity;
            let closest_obj = null;

            for (const mesh of scene.meshes) {
                let t = intersectRayTriangle(origin, direction, mesh);
                if (t && t >= min_t && t <= max_t && t < closest_t) {
                    closest_t = t;
                    closest_obj = mesh;
                }
            }
            
            for (const sphere of scene.spheres) {
                let t = intersectRaySphere(origin, direction, sphere);
                if (t && t >= min_t && t <= max_t && t < closest_t) {
                    closest_t = t;
                    closest_obj = sphere;
                }
            }

            return closest_obj ? { obj: closest_obj, t: closest_t } : null;
        }

        function computeLighting(point, normal, view, specular, scene) {
            let intensity = 0;

            for (const light of scene.lights) {
                if (light.type === ${LightType.AMBIENT}) {
                    intensity += light.intensity;
                    continue;
                }

                let light_vec, t_max;
                if (light.type === ${LightType.POINT}) {
                    light_vec = light.position.sub(point);
                    t_max = 1;
                } else { // DIRECTIONAL
                    light_vec = light.position;
                    t_max = Infinity;
                }

                // Shadow check
                const blocker = closestIntersection(point, light_vec, EPSILON, t_max, scene);
                if (blocker) continue;

                // Diffuse
                const n_dot_l = normal.dot(light_vec);
                if (n_dot_l > 0) {
                    intensity += light.intensity * n_dot_l / (normal.length() * light_vec.length());
                }

                // Specular
                if (specular > 0) {
                    const r = ReflectRay(light_vec, normal);
                    const r_dot_v = r.dot(view);
                    if (r_dot_v > 0) {
                        intensity += light.intensity * Math.pow(r_dot_v / (r.length() * view.length()), specular);
                    }
                }
            }

            return Math.min(1, intensity);
        }

        onmessage = function (e) {
            const data = e.data;
            const sceneData = data.scene;
            
            // Reconstruct scene with class instances
            const scene = {
                ...sceneData,
                camera_position: new Vec(sceneData.camera_position.x, sceneData.camera_position.y, sceneData.camera_position.z),
                background_color: new Color(sceneData.background_color.r, sceneData.background_color.g, sceneData.background_color.b),
                meshes: sceneData.meshes.map(mesh => {
                    const vertices = mesh.vertices.map(v => new Vec(v.x, v.y, v.z));
                    const triangle = {
                        type: 'triangle',
                        vertices,
                        normal: mesh.normal ? new Vec(mesh.normal.x, mesh.normal.y, mesh.normal.z) : null,
                        color: new Color(mesh.color.r, mesh.color.g, mesh.color.b),
                        material: {
                            diffuse: new Color(
                                mesh.material.diffuse.r,
                                mesh.material.diffuse.g,
                                mesh.material.diffuse.b
                            ),
                            specular: mesh.material.specular,
                            reflective: mesh.material.reflective
                        }
                    };
                    if (!triangle.normal) {
                        const edge1 = vertices[1].sub(vertices[0]);
                        const edge2 = vertices[2].sub(vertices[0]);
                        triangle.normal = edge1.cross(edge2).normalized();
                    }
                    return triangle;
                }),
                spheres: sceneData.spheres.map(s => ({
                    type: s.type,
                    center: new Vec(s.center.x, s.center.y, s.center.z),
                    radius: s.radius,
                    color: new Color(s.color.r, s.color.g, s.color.b),
                    specular: s.specular,
                    reflective: s.reflective,
                    material: {
                        diffuse: new Color(s.color.r, s.color.g, s.color.b),
                        specular: s.specular,
                        reflective: s.reflective
                    }
                })),
                lights: sceneData.lights.map(light => ({
                    ...light,
                    position: light.position ? new Vec(light.position.x, light.position.y, light.position.z) : null
                }))
            };

            const halfWidth = Math.floor(data.width / 2);
            const pixels = [];

            for (let y = data.startY; y < data.endY; y++) {
                for (let x = -halfWidth; x < halfWidth; x++) {
                    const direction = canvasToViewport(x, y, data.width, data.height, scene);
                    const color = traceRay(
                        scene.camera_position,
                        direction,
                        0.5,
                        Infinity,
                        scene.recursion_limit,
                        scene
                    );
                    pixels.push({ x, y, color });
                }
            }

            postMessage({ pixels });
        };
    `;

    // Create a blob URL for the worker
    const blob = new Blob([workerCode], { type: 'application/javascript' });
    const workerUrl = URL.createObjectURL(blob);

    for (let i = 0; i < numWorkers; i++) {
        const worker = new Worker(workerUrl);
        const startY = -Math.floor(height / 2) + i * rowsPerWorker;
        const endY = Math.min(startY + rowsPerWorker, Math.floor(height / 2));

        worker.postMessage({
            scene: {
                // Serialize scene data
                viewport_size: scene.viewport_size,
                projection_plane_z: scene.projection_plane_z,
                camera_position: {x: scene.camera_position.x, y: scene.camera_position.y, z: scene.camera_position.z},
                background_color: {r: scene.background_color.r, g: scene.background_color.g, b: scene.background_color.b},
                recursion_limit: scene.recursion_limit,
                meshes: scene.meshes.map(mesh => ({
                    type: mesh.type,
                    vertices: mesh.vertices.map(v => ({x: v.x, y: v.y, z: v.z})),
                    normal: mesh.normal ? {x: mesh.normal.x, y: mesh.normal.y, z: mesh.normal.z} : null,
                    color: {r: mesh.color.r, g: mesh.color.g, b: mesh.color.b},
                    material: {
                        diffuse: {r: mesh.material.diffuse.r, g: mesh.material.diffuse.g, b: mesh.material.diffuse.b},
                        specular: mesh.material.specular,
                        reflective: mesh.material.reflective
                    }
                })),
                spheres: scene.spheres.map(s => ({
                    type: s.type,
                    center: new Vec(s.center.x, s.center.y, s.center.z),
                    radius: s.radius,
                    color: new Color(s.color.r, s.color.g, s.color.b),
                    specular: s.specular,
                    reflective: s.reflective,
                    material: {
                        diffuse: new Color(s.color.r, s.color.g, s.color.b),
                        specular: s.specular,
                        reflective: s.reflective
                    }
                })),
                lights: scene.lights.map(light => ({
                    type: light.type,
                    intensity: light.intensity,
                    position: light.position ? {x: light.position.x, y: light.position.y, z: light.position.z} : null
                }))
            },
            width,
            height,
            startY,
            endY
        });

        worker.onmessage = function (e) {
            const { pixels } = e.data;
            for (const { x, y, color } of pixels) {
                putPixel(x, y, color);
            }
            completedWorkers++;
            if (completedWorkers === numWorkers) {
                updateCanvas();
                // Clean up worker URLs
                URL.revokeObjectURL(workerUrl);
            }
        };
    }
}

function renderSpheres() {
    scene.meshes = []; // Clear meshes if rendering spheres
    scene.spheres = [
        new Sphere(new Vec(0, -1, 3), 1, new Color(255, 0, 0), 500, 0.2),
        new Sphere(new Vec(-2, 0, 4), 1, new Color(0, 255, 0), 10, 0.4),
        new Sphere(new Vec(2, 0, 4), 1, new Color(0, 0, 255), 500, 0.3),
        new Sphere(new Vec(0, -5001, 0), 5000, new Color(255, 255, 0), 1000, 0.5)
    ];
    render();
}

// Initial render
render();