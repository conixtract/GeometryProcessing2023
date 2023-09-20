#include "MeshRenderable.hh"

#include <glow/common/scoped_gl.hh>
#include <glow/objects/TextureCubeMap.hh>

#include "../RenderInfo.hh"
#include "../detail/MeshShaderBuilder.hh"

using namespace glow;
using namespace glow::viewer;

aabb MeshRenderable::computeAabb() { return getMeshAABB().transformed(transform()); }

void MeshRenderable::renderMesh(const RenderInfo& info)
{
    if (mMesh->getVertexCount() == 0)
        return; // empty

    auto shader = mForwardShader->use();

    GLOW_SCOPED(onlyEnableIf, getBackfaceCullingEnabled(), GL_CULL_FACE);

    shader["uModel"] = transform();
    shader["uNormalModel"] = inverse(transpose(tg::mat3(transform())));
    shader["uView"] = info.view;
    shader["uProj"] = info.proj;
    shader["uFresnel"] = getFresnel();
    shader["uIsTransparent"] = getRenderMode() == RenderMode::Transparent;
    shader["uIsShadingEnabled"] = getShadingEnabled();
    shader["uSeed"] = tg::u32(info.accumulationCount);
    shader["uCamPos"] = info.camPos;
    shader["uClipPlane"] = getClipPlane();

    if (getEnvMap())
    {
        shader["uEnvMap"] = getEnvMap();
        shader["uEnvRefl"] = getEnvReflectivity();
    }

    if (getColorMapping())
        getColorMapping()->prepareShader(shader);
    if (getTexturing())
        getTexturing()->prepareShader(shader);
    if (getMasking())
        getMasking()->prepareShader(shader);
    for (auto const& a : getAttributes())
        a->prepareShader(shader);
    mMesh->bind().draw();
}

void MeshRenderable::renderShadow(RenderInfo const& info) { renderMesh(info); }

void MeshRenderable::renderForward(RenderInfo const& info)
{
    if (getRenderMode() != RenderMode::Opaque)
        return;

    renderMesh(info);
}

void MeshRenderable::renderTransparent(RenderInfo const& info)
{
    if (getRenderMode() != RenderMode::Transparent)
        return;

    renderMesh(info);
}

/*picking functionality*/
void MeshRenderable::renderPicking(RenderInfo const& info, int32_t renderableID)
{
    if (!hasPicker() || (mMesh->getVertexCount() == 0))
        return;

    auto shader = mPickingShader->use();

    GLOW_SCOPED(onlyEnableIf, getBackfaceCullingEnabled(), GL_CULL_FACE);

    shader["uModel"] = transform();
    shader["uNormalModel"] = inverse(transpose(tg::mat3(transform())));
    shader["uView"] = info.view;
    shader["uProj"] = info.proj;
    shader["uRenderableID"] = renderableID;

    mMesh->bind().draw();
}

SharedMeshRenderable MeshRenderable::create(builder::PolygonBuilder const& builder)
{
    auto r = std::make_shared<MeshRenderable>();
    r->initGeometry(builder.getMeshDef(), builder.getAttributes());
    return r;
}

void MeshRenderable::init()
{
    // add missing attributes
    if (!hasAttribute("aNormal"))
        addAttribute(getMeshDefinition()->computeFaceNormalAttribute());
    if (!hasAttribute("aColor"))
        addAttribute(detail::make_mesh_attribute("aColor", tg::color3::white));
    if (getMasking())
        addAttribute(getMasking()->mDataAttribute);
    if (getTexturing())
        addAttribute(getTexturing()->mCoordsAttribute);
    else if (getColorMapping())
        addAttribute(getColorMapping()->mDataAttribute);

    if (hasPicker() && !hasAttribute("aPickID"))
    {
        /* Define ids per primitive*/
        pm::face_attribute<int32_t> IDs;
        int32_t vID = 0;
        auto& p = this->getPicker();

        auto pm = std::dynamic_pointer_cast<detail::PolyMeshDefinition>(getMeshDefinition());
        TG_ASSERT(pm);

        IDs = pm::face_attribute<int32_t>(pm->mesh);
        std::vector<pm::face_index> face_indices;
        face_indices.reserve(IDs.size());
        for (auto f : pm->mesh.faces())
        {
            IDs[f] = vID;
            vID++;
            face_indices.push_back(f.idx); // offers ability to map back the read IDs to face_indices
        }

        addAttribute(detail::make_mesh_attribute("aPickID", IDs));
        p.init(face_indices);
    }

    const auto aColor = getAttribute("aColor");

    // build meshes
    {
        std::vector<SharedArrayBuffer> abs;
        SharedArrayBuffer posAB;
        SharedElementArrayBuffer eab;

        for (auto const& attr : getAttributes())
        {
            auto ab = attr->createMeshRenderableArrayBuffer(*getMeshDefinition());
            if (!ab)
                continue;


            abs.push_back(ab);
        }

        mMesh = VertexArray::create(abs, getIndexBuffer());
    }

    // build shader
    {
        detail::MeshShaderBuilder sb;

        sb.addUniform("mat4", "uModel");
        sb.addUniform("mat4", "uProj");
        sb.addUniform("mat4", "uView");
        sb.addUniform("mat3", "uNormalModel");
        sb.addUniform("vec3", "uCamPos");
        sb.addUniform("bool", "uIsTransparent");
        sb.addUniform("bool", "uFresnel");
        sb.addUniform("bool", "uIsShadingEnabled");
        sb.addUniform("uint", "uSeed");
        sb.addUniform("vec4", "uClipPlane");

        if (getEnvMap())
        {
            sb.addUniform("samplerCube", "uEnvMap");
            sb.addUniform("float", "uEnvRefl");
        }

        sb.addFragmentLocation("vec4", "fColor");
        sb.addFragmentLocation("vec3", "fNormal");

        sb.addPassthrough("vec3", "Normal");
        sb.addPassthrough("vec3", "WorldPos");
        sb.addPassthrough("float", "VertexID");
        sb.addPassthrough(aColor->typeInShader(), "Color", detail::MeshShaderBuilder::TypeHandling::ExtendToVec4Color);

        for (auto const& attr : getAttributes())
            attr->buildShader(sb);

        // data mapped mesh
        if (getColorMapping())
            getColorMapping()->buildShader(sb);

        // texture mesh
        if (getTexturing())
            getTexturing()->buildShader(sb);

        // masked mesh
        if (getMasking())
            getMasking()->buildShader(sb);

        sb.addVertexShaderCode(R"(
                                gl_Position = uProj * uView * uModel * vec4(aPosition, 1.0);
                                vOut.VertexID = float(gl_VertexID);
                                vOut.WorldPos = vec3(uModel * vec4(aPosition, 1.0));
                               )");

        sb.addFragmentShaderCode(R"(
                                    if (dot(uClipPlane.xyz, vWorldPos) > uClipPlane.w)
                                        discard;

                                    fNormal = normalize(uNormalModel * vNormal) * (gl_FrontFacing ? 1 : -1);
                                    if(uIsShadingEnabled)
                                    {
                                        fColor.rgb = vColor.rgb * (fNormal.y * .4 + .6);
                                    }
                                    else
                                    {
                                        fColor.rgb = vColor.rgb;
                                    }
                                    fColor.a = 1;
                                 )");

        if (getEnvMap())
        {
            sb.addFragmentShaderCode(R"(
                                     vec3 V = normalize(uCamPos - vIn.WorldPos);
                                     vec3 R = reflect(-V, fNormal);
                                     fColor.rgb += texture(uEnvMap, R).rgb * uEnvRefl;
                                     )");
        }

        sb.addFragmentShaderCode(R"(
                                    fColor = clamp(fColor, 0., 1.);
                                 )");


        if (getRenderMode() == RenderMode::Transparent)
            sb.addFragmentShaderCode(R"(
                                     if (uIsTransparent)
                                     {
                                        uint h = uint(gl_FragCoord.x) * 4096 + uint(gl_FragCoord.y);
                                        h = hash_combine(h, uint(uSeed));
                                        h = hash_combine(h, uint(vVertexID * 17));
                                        h = wang_hash(h);

                                        float a = vColor.a;

                                        if (uFresnel)
                                        {
                                            vec3 V = normalize(uCamPos - vWorldPos);
                                            float t = 1 - abs(dot(fNormal, V));
                                            t = (t * t) * (t * t) * t;
                                            a = mix(vColor.a, 1, t);
                                        }

                                        if (wang_float(h) > a)
                                            discard;
                                     })");

        mForwardShader = sb.createProgram();

        // create Picking Shader -> use
        if (hasPicker())
        {
            detail::MeshShaderBuilder sb;

            sb.addUniform("mat4", "uModel");
            sb.addUniform("mat4", "uProj");
            sb.addUniform("mat4", "uView");
            sb.addUniform("int", "uRenderableID");

            sb.addFragmentLocation("ivec2", "fPickIDs");

            sb.addPassthrough("int", "FragmentID");
            sb.addPassthrough("int", "RenderableID");

            for (auto const& attr : getAttributes())
                attr->buildShader(sb);

            sb.addVertexShaderCode(R"(
                                gl_Position = uProj * uView * uModel * vec4(aPosition, 1.0);
                                vOut.FragmentID = aPickID;
                                vOut.RenderableID = uRenderableID;
                               )");

            sb.addFragmentShaderCode(R"(
                                  fPickIDs = ivec2((vIn.RenderableID), (vIn.FragmentID));
                                 )");

            mPickingShader = sb.createProgram();
        }
    }
}

size_t MeshRenderable::computeHash() const { return computeGenericGeometryHash(); }
