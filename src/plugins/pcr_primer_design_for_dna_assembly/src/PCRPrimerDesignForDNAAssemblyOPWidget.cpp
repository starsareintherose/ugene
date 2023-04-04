/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2021 UniPro <ugene@unipro.ru>
 * http://ugene.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include "PCRPrimerDesignForDNAAssemblyOPWidget.h"
#include "tasks/ExtractPrimerTask.h"
#include "tasks/PCRPrimerDesignForDNAAssemblyTask.h"

#include <U2Algorithm/TmCalculatorRegistry.h>

#include <U2Core/AppContext.h>
#include <U2Core/AppSettings.h>
#include <U2Core/BaseDocumentFormats.h>
#include <U2Core/CreateAnnotationTask.h>
#include <U2Core/DNAAlphabet.h>
#include <U2Core/DNASequenceObject.h>
#include <U2Core/DNASequenceSelection.h>
#include <U2Core/DNASequenceUtils.h>
#include <U2Core/DocumentModel.h>
#include <U2Core/FileFilters.h>
#include <U2Core/GObjectTypes.h>
#include <U2Core/GUrlUtils.h>
#include <U2Core/IOAdapter.h>
#include <U2Core/L10n.h>
#include <U2Core/PrimerValidator.h>
#include <U2Core/U2DbiRegistry.h>
#include <U2Core/U2OpStatusUtils.h>
#include <U2Core/ProjectModel.h>
#include <U2Core/U2SafePoints.h>
#include <U2Core/UserApplicationsSettings.h>

#include <U2Gui/CreateAnnotationWidgetController.h>
#include <U2Gui/LastUsedDirHelper.h>
#include <U2Gui/PrimerLineEdit.h>
#include <U2Gui/U2FileDialog.h>
#include <U2Gui/U2WidgetStateStorage.h>

#include <U2View/ADVSequenceObjectContext.h>
#include <U2View/AnnotatedDNAView.h>

#include <QVBoxLayout>

namespace {
// When a non-nucleotide sequence is selected, the widget should be disabled. But the groupbox titles continue to be
// black as if they enabled. This code makes them gray.
//void makeGroupboxTittleGrayIfDisable(QGroupBox *gb) {
//    QPalette palette;
//    palette.setColor(QPalette::Disabled, QPalette::WindowText, QApplication::palette().color(QPalette::Disabled,
//                                                                                             QPalette::WindowText));
//    gb->setPalette(palette);
//}
}  // namespace

namespace U2 {

const QString PCRPrimerDesignForDNAAssemblyOPWidget::USER_PRIMERS_SHOW_HIDE_ID = "user-primers-show-hide-id";
const QString PCRPrimerDesignForDNAAssemblyOPWidget::PICK_PRIMERS_SETTINGS_HIDE_ID = "pick-primers-settings-show-hide-id";
const QString PCRPrimerDesignForDNAAssemblyOPWidget::TM_SETTINGS_SHOW_HIDE_ID = "tm-settings-show-hide-id";
const QString PCRPrimerDesignForDNAAssemblyOPWidget::PCR_TABLE_OBJECT_NAME =
    QApplication::translate("PCRPrimerDesignForDNAAssemblyOPWidget", "PCR Primer Design for DNA assembly");

static constexpr char* ID_POSTFIX = "_pcr_primer_design";

PCRPrimerDesignForDNAAssemblyOPWidget::PCRPrimerDesignForDNAAssemblyOPWidget(AnnotatedDNAView* _annDnaView)
    : QWidget(nullptr),
      annDnaView(_annDnaView),
      savableWidget(this, GObjectViewUtils::findViewByName(annDnaView->getName())),
      tmCalculatorId(annDnaView->getName() + ID_POSTFIX) {
    setupUi(this);
    parametersMinMaxSpinBoxes = { { sbMinRequireGibbs, sbMaxRequireGibbs },
                                  { spMinRequireMeltingTeml, spMaxRequireMeltingTeml },
                                  { spMinRequirePrimerLength, spMaxRequirePrimerLength } };

    userPrimersShowHideGroup->init(USER_PRIMERS_SHOW_HIDE_ID, tr("Check primers"), wgtUserPrimers, false);
    mainSettingsShowHide->init(PICK_PRIMERS_SETTINGS_HIDE_ID, tr("Pick primers"), wgtMainSettings, true);
    {
        auto tempCalculator = AppContext::getTmCalculatorRegistry()->createTmCalculator(tmCalculatorId);
        wgtTmSettings->init(tempCalculator->getSettings());
    }
    wgtTmSettingsShowHide->init(TM_SETTINGS_SHOW_HIDE_ID, tr("Melting temperature"), wgtTmSettings, false);

    {
        CreateAnnotationModel createAnnotationModel;
        createAnnotationModel.data->name = "top_primers";
        createAnnotationModel.sequenceObjectRef = annDnaView->getActiveSequenceContext()->getSequenceObject()->getReference();
        createAnnotationModel.hideAnnotationType = true;
        createAnnotationModel.hideAnnotationName = false;
        createAnnotationModel.hideLocation = true;
        annWgtController = new CreateAnnotationWidgetController(createAnnotationModel, this, CreateAnnotationWidgetController::OptionsPanel);
        wgtShowHideMain->layout()->addWidget(annWgtController->getWidget());
    }


    auto seqLength = annDnaView->getActiveSequenceContext()->getSequenceLength();

    // Default values, have been provided by the customer
    // Left area: from "1" to "10% of sequence length"
    // Right area: from "50% of sequence length" to "60% of sequence length"
    sbLeftAreaStart->setValue(1);
    sbLeftAreaEnd->setValue(seqLength / 10);
    sbRightAreaStart->setValue(seqLength / 10 * 5);
    sbRightAreaEnd->setValue(seqLength / 10 * 6);

    //makeGroupboxTittleGrayIfDisable(groupBox);
    //makeGroupboxTittleGrayIfDisable(groupBox_2);
    makeWarningInvisibleIfDna();

    connect(pbStart, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_start);
    connect(tbLeftAreaSelectManually, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_selectManually);
    connect(tbRightAreaSelectManually, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_selectManually);
    U2WidgetStateStorage::restoreWidgetState(savableWidget);
    initAreaSpinboxes();  // Must be after sbArea->setValue and restoreWidgetState.

    sl_updateParametersRanges();
    const auto& parametersMinMaxSpinBoxesKeys = parametersMinMaxSpinBoxes.keys();
    for (auto minSb : qAsConst(parametersMinMaxSpinBoxesKeys)) {
        connect(minSb, &QSpinBox::editingFinished, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_minValueChanged);
        auto maxSb = parametersMinMaxSpinBoxes.value(minSb);
        SAFE_POINT(maxSb != nullptr, L10N::nullPointerError("QSpinBox"), );

        connect(maxSb, &QSpinBox::editingFinished, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_maxValueChanged);
    }
    sl_addEndingOnOverhangLabels(cbOverhangConnection->currentIndex());
    connect(cbOverhangConnection, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_addEndingOnOverhangLabels);
    connect(tbLoadOtherSequencesInPcr, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_loadOtherSequenceInPcr);
    connect(annDnaView, SIGNAL(si_activeSequenceWidgetChanged(ADVSequenceWidget *, ADVSequenceWidget *)),
            SLOT(sl_activeSequenceChanged()));
    connect(annDnaView, &AnnotatedDNAView::si_sequenceModified, this,
            &PCRPrimerDesignForDNAAssemblyOPWidget::sl_sequenceModified);
    connect(annDnaView->getActiveSequenceContext()->getSequenceObject(), &U2SequenceObject::si_sequenceChanged, this,
            &PCRPrimerDesignForDNAAssemblyOPWidget::sl_sequenceModified);

    connect(pbFindReverseComplement, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_selectReverseComplementInTable);
    connect(pbAddToForward5, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_add5ForwardSequence);
    connect(pbAddToForward3, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_add3ForwardSequence);
    connect(pbAddToReverse5, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_add5ReverseSequence);
    connect(pbAddToReverse3, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_add3ReverseSequence);

    leFilter->setValidator(new PrimerValidator(this, false));
    connect(leFilter, &QLineEdit::textEdited, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_updateSequenceList);
    connect(pbAddToReverse3, &QAbstractButton::clicked, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_add3ReverseSequence);
}

PCRPrimerDesignForDNAAssemblyOPWidget::~PCRPrimerDesignForDNAAssemblyOPWidget() {
    AppContext::getTmCalculatorRegistry()->saveSettings(tmCalculatorId, wgtTmSettings->getSettings());
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_activeSequenceChanged() {
    makeWarningInvisibleIfDna();
    sl_sequenceModified();

    // Release "Select manually" buttons.
    tbLeftAreaSelectManually->setChecked(false);
    tbRightAreaSelectManually->setChecked(false);
    smButton = nullptr;
    sbStartRegion = nullptr;
    sbEndRegion = nullptr;
    disconnect(updateRegionConnection);
    updateRegionConnection = QMetaObject::Connection();
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_sequenceModified() {
    ADVSequenceObjectContext *sequenceContext = annDnaView->getActiveSequenceContext();
    CHECK(sequenceContext != nullptr, )
    qint64 seqLength = sequenceContext->getSequenceLength();

    // End spinbox maximum value is the previous sequence length. If it has changed, need to update spinboxes.
    if (seqLength != sbLeftAreaEnd->maximum()) {
        int start = sbLeftAreaStart->value() - 1;
        int len = sbLeftAreaEnd->value() - start;
        setRegion(sbLeftAreaStart, sbLeftAreaEnd, {start, len}, seqLength);

        start = sbRightAreaStart->value() - 1;
        len = sbRightAreaEnd->value() - start;
        setRegion(sbRightAreaStart, sbRightAreaEnd, {start, len}, seqLength);

        sl_updateParametersRanges();
    }
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_start() {
    PCRPrimerDesignForDNAAssemblyTaskSettings settings;
    settings.forwardUserPrimer = leForwardPrimer->text();
    settings.reverseUserPrimer = leReversePrimer->text();

    settings.gibbsFreeEnergy.minValue = sbMinRequireGibbs->value();
    settings.gibbsFreeEnergy.maxValue = sbMaxRequireGibbs->value();
    settings.meltingPoint.minValue = spMinRequireMeltingTeml->value();
    settings.meltingPoint.maxValue = spMaxRequireMeltingTeml->value();
    settings.primerLength.minValue = spMinRequirePrimerLength->value();
    settings.primerLength.maxValue = spMaxRequirePrimerLength->value();

    settings.minGibbs = sbMinGibbs->value();
    settings.maxTm = spMaxTm->value();
    settings.maxLength = spMaxLength->value();

    settings.insertTo = static_cast<PCRPrimerDesignForDNAAssemblyTaskSettings::OverhangConnection>(cbOverhangConnection->currentIndex());
    settings.leftPrimerOverhang = leLeftOverhang->text();
    settings.rightPrimerOverhang = leRightOverhang->text();

    settings.leftArea.startPos = sbLeftAreaStart->value() - 1;
    settings.leftArea.length = sbLeftAreaEnd->value() - sbLeftAreaStart->value();
    settings.rightArea.startPos = sbRightAreaStart->value() - 1;
    settings.rightArea.length = sbRightAreaEnd->value() - sbRightAreaStart->value();

    settings.otherSequencesInPcrUrl = leOtherSequencesInPcrFilePath->text();

    settings.searchForAdditionslPrimers = cbAdditionslPrimers->isChecked();

    auto settingsMap = wgtTmSettings->getSettings();
    auto id = settingsMap.value(TmCalculator::KEY_ID).toString();
    settings.tmCalculator = AppContext::getTmCalculatorRegistry()->getById(id)->createCalculator(settingsMap);

    auto activeSequenceContext = annDnaView->getActiveSequenceContext();
    SAFE_POINT(activeSequenceContext != nullptr, L10N::nullPointerError("ADVSequenceObjectContext"), );

    U2SequenceObject* sequenceObject = activeSequenceContext->getSequenceObject();
    U2OpStatus2Log os;
    auto sequence = sequenceObject->getWholeSequenceData(os);
    CHECK_OP(os, );

    qint64 sequenceLength = activeSequenceContext->getSequenceLength();
    auto checkArea = [sequenceLength](const U2Region &area) {
        return area.startPos >= 0 && area.length >= 0 && area.endPos() < sequenceLength;
    };
    SAFE_POINT(checkArea(settings.leftArea) && checkArea(settings.rightArea),
               "Invalid areas, PCR Primer Design For DNA Assembly Task cannot be started", )

    pcrTask = new PCRPrimerDesignForDNAAssemblyTask(settings, sequence);
    connect(pcrTask, SIGNAL(si_stateChanged()), SLOT(sl_onFindTaskFinished()));
    AppContext::getTaskScheduler()->registerTopLevelTask(pcrTask);
    pbStart->setEnabled(false);
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_selectManually() {
    auto currentSMButton = qobject_cast<QToolButton*>(sender());
    SAFE_POINT(currentSMButton != nullptr, L10N::nullPointerError("QAbstractButton"), );

    if (smButton != nullptr) {
        bool disableState = true;
        if (smButton->isChecked()) {
            smButton->setChecked(false);
            disableState = false;
        }
        smButton = nullptr;
        sbStartRegion = nullptr;
        sbEndRegion = nullptr;
        disconnect(updateRegionConnection);
        updateRegionConnection = QMetaObject::Connection();

        if (disableState) {
            return;
        }
    }

    smButton = currentSMButton;
    if (smButton == tbLeftAreaSelectManually) {
        sbStartRegion = sbLeftAreaStart;
        sbEndRegion = sbLeftAreaEnd;
    } else {
        sbStartRegion = sbRightAreaStart;
        sbEndRegion = sbRightAreaEnd;
    }
    auto sequenceSelection = annDnaView->getActiveSequenceContext()->getSequenceSelection();
    int startRegionValue = sbStartRegion->value() - 1;
    int endRegionValue = sbEndRegion->value();
    sequenceSelection->setSelectedRegions({U2Region(startRegionValue, endRegionValue - startRegionValue)});
    updateRegionConnection = connect(sequenceSelection, &DNASequenceSelection::si_onSelectionChanged, this, &PCRPrimerDesignForDNAAssemblyOPWidget::sl_updateRegion);
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_updateRegion() {
    const auto& selectedRegions = qobject_cast<DNASequenceSelection*>(sender())->getSelectedRegions();
    CHECK(!selectedRegions.isEmpty(), );

    if (selectedRegions.size() > 1) {
        coreLog.details(tr("Multiple regions selection, the only first is set as primer search area."));
    }
    const auto& region = selectedRegions.first();
    SAFE_POINT(sbStartRegion != nullptr, L10N::nullPointerError("QSpinBox"), );

    setRegion(sbStartRegion, sbEndRegion, region);
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_updateParametersRanges() {
    auto seqLength = annDnaView->getActiveSequenceContext()->getSequenceLength();
    auto parametersMinMaxSpinBoxesKeys = parametersMinMaxSpinBoxes.keys();
    for (const auto& minSb : qAsConst(parametersMinMaxSpinBoxesKeys)) {
        auto maxSb = parametersMinMaxSpinBoxes.value(minSb);
        SAFE_POINT(maxSb != nullptr, L10N::nullPointerError("QSpinBox"), );

        minSb->setMaximum(maxSb->value() - 1);
        maxSb->setMaximum(seqLength);
        maxSb->setMinimum(minSb->value() + 1);
    }
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_minValueChanged() {
    auto minSb = qobject_cast<QSpinBox*>(sender());
    SAFE_POINT(minSb != nullptr, L10N::nullPointerError("QSpinBox"), );

    auto maxSb = parametersMinMaxSpinBoxes.value(minSb);
    SAFE_POINT(maxSb != nullptr, L10N::nullPointerError("QSpinBox"), );

    maxSb->setMinimum(minSb->value() + 1);
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_maxValueChanged() {
    auto maxSb = qobject_cast<QSpinBox*>(sender());
    SAFE_POINT(maxSb != nullptr, L10N::nullPointerError("QSpinBox"), );

    auto minSb = parametersMinMaxSpinBoxes.key(maxSb);
    SAFE_POINT(minSb != nullptr, L10N::nullPointerError("QSpinBox"), );

    minSb->setMaximum(maxSb->value() - 1);
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_addEndingOnOverhangLabels(int index) {
    static const QString DIRECT = QT_TR_NOOP("Info: this overhang is located on direct strand");
    static const QString REV_COMP = QT_TR_NOOP("Info: this overhang is located on reverse-complementary strand (opposite direction)");

    if (index == 0) {
        lbLeftOverhangInfo->setText(DIRECT);
        lbRightOverhangInfo->setText(REV_COMP);
    } else {
        lbLeftOverhangInfo->setText(REV_COMP);
        lbRightOverhangInfo->setText(DIRECT);
    }
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_loadOtherSequenceInPcr() {
    QString filter = FileFilters::createFileFilterByObjectTypes({ GObjectTypes::SEQUENCE }, true);
    QString selectedFilter = FileFilters::createFileFilterByObjectTypes({ GObjectTypes::SEQUENCE });
    LastUsedDirHelper lod;
    QString file = U2FileDialog::getOpenFileName(nullptr, tr("Select an \"Other sequences in PCR reaction\" file"), lod, filter, selectedFilter);
    CHECK(!file.isEmpty(), );

    leOtherSequencesInPcrFilePath->setText(file);

}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_onFindTaskFinished() {
    CHECK(sender() == pcrTask, );
    SAFE_POINT(nullptr != pcrTask, L10N::nullPointerError("InSilicoPcrTask"), );
    if (pcrTask->isCanceled() || pcrTask->hasError()) {
        pbStart->setEnabled(true);
        disconnect(pcrTask, SIGNAL(si_stateChanged()));
        pcrTask = nullptr;
        return;
    }
    CHECK(pcrTask->isFinished(), );
    pbStart->setEnabled(true);
    productsTable->setCurrentProducts(pcrTask->getResults(), annDnaView);
    createResultAnnotations();
    backboneSequence = pcrTask->getBackboneSequence();
    lastRunSettings = pcrTask->getSettings();
    pcrTask = nullptr;
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_selectReverseComplementInTable() {
    auto sequence = getSelectedSequence();
    auto dnaAlphabet = AppContext::getDNAAlphabetRegistry()->findById(BaseDNAAlphabetIds::NUCL_DNA_DEFAULT());
    SAFE_POINT(dnaAlphabet != nullptr, L10N::nullPointerError("DNAAlphabet"), );

    auto reverseComplementSequence = DNASequenceUtils::reverseComplement(sequence.toLocal8Bit(), dnaAlphabet);
    for (int i = 0; i < twGeneratedSequences->rowCount(); i++) {
        auto txt = twGeneratedSequences->item(i, 0)->text();
        CHECK_CONTINUE(txt == reverseComplementSequence);

        twGeneratedSequences->selectRow(i);
        twGeneratedSequences->setFocus();
        break;
    }
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_add5ForwardSequence() {
    leForwardPrimer->setLeftEnd(getSelectedSequence());
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_add3ForwardSequence() {
    leForwardPrimer->setRightEnd(getSelectedSequence());

}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_add5ReverseSequence() {
    leReversePrimer->setRightEnd(getSelectedSequence());
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_add3ReverseSequence() {
    leReversePrimer->setLeftEnd(getSelectedSequence());
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_updateSequenceList(const QString& text) {
    twGeneratedSequences->updateSequenceList(text);
}

void PCRPrimerDesignForDNAAssemblyOPWidget::sl_annotationCreationTaskFinished() {
    auto t = qobject_cast<CreateAnnotationsTask *>(sender());
    if (t->getState() != Task::State_Finished || t->isCanceled() || t->hasError()) {
        return;
    }
    CHECK(t->getResultAnnotations().size() > 0, );
    productsTable->setAnnotationGroup(t->getResultAnnotations().at(0)->getGroup());
}

void PCRPrimerDesignForDNAAssemblyOPWidget::createResultAnnotations() {
    annWgtController->prepareAnnotationObject();
    const CreateAnnotationModel& model = annWgtController->getModel();
    auto ato = model.getAnnotationObject();

    QMap<QString, QList<SharedAnnotationData>> resultAnnotations;
    auto results = pcrTask->getResults();
    for (int i = 0; i < results.size(); i++) {
        const auto& regions = results[i];
        CHECK_CONTINUE(!regions.first.isEmpty());

        QList<SharedAnnotationData> annotations;
        auto createAnnotation = [](const U2Region& region, U2Strand strand) -> SharedAnnotationData {
            SharedAnnotationData data(new AnnotationData());
            data->setStrand(strand);
            data->name = "pcr_design_primers";
            data->type = U2FeatureTypes::Primer;
            data->location->regions.append(region);
            //U2Qualifier les("left_end_seq",
            //data->qualifiers
            return data;
        };
        annotations << createAnnotation(regions.first, U2Strand::Direct);
        annotations << createAnnotation(regions.second, U2Strand::Complementary);
        resultAnnotations["pcr_design_primers/pair " + QString::number(i + 1)].append(annotations);

    }
    CHECK(!resultAnnotations.isEmpty(), );

    AppContext::getTaskScheduler()->registerTopLevelTask(new CreateAnnotationsTask(ato, resultAnnotations));
}

void PCRPrimerDesignForDNAAssemblyOPWidget::makeWarningInvisibleIfDna() {
    ADVSequenceObjectContext *sequenceContext = annDnaView->getActiveSequenceContext();
    CHECK(sequenceContext != nullptr, )

    bool isDna = sequenceContext->getAlphabet()->isDNA();
    runPcrPrimerDesignWidget->setEnabled(isDna);
    alphabetWarningLabel->setVisible(!isDna);
}

void PCRPrimerDesignForDNAAssemblyOPWidget::setRegion(QSpinBox *start, QSpinBox* end, const U2Region &region) const {
    ADVSequenceObjectContext *sequenceContext = annDnaView->getActiveSequenceContext();
    SAFE_POINT(start != nullptr && end != nullptr, L10N::nullPointerError("QSpinBox"), )
    CHECK(sequenceContext != nullptr, )

    qint64 seqLength = sequenceContext->getSequenceLength();
    setRegion(start, end, region, seqLength);
}

void PCRPrimerDesignForDNAAssemblyOPWidget::setRegion(QSpinBox *start,
                                                      QSpinBox *end,
                                                      const U2Region &region,
                                                      qint64 seqLength) const {
    start->setRange(1, seqLength);
    end->setRange(1, seqLength);

    start->setValue(region.startPos + 1);
    end->setValue(region.endPos());
    updateAreaSpinboxRanges(start, end, seqLength);
}

void PCRPrimerDesignForDNAAssemblyOPWidget::updateAreaSpinboxRanges(QSpinBox *start,
                                                                    QSpinBox *end,
                                                                    qint64 seqLength) const {
    end->setMaximum(seqLength);
    start->setMaximum(end->value());
    end->setMinimum(start->value());
}

void PCRPrimerDesignForDNAAssemblyOPWidget::initAreaSpinboxes() const {
    ADVSequenceObjectContext *context = annDnaView->getActiveSequenceContext();
    SAFE_POINT(context != nullptr, L10N::nullPointerError("ADVSequenceObjectContext"), )
    qint64 sequenceLength = context->getSequenceLength();
    updateAreaSpinboxRanges(sbLeftAreaStart, sbLeftAreaEnd, sequenceLength);
    updateAreaSpinboxRanges(sbRightAreaStart, sbRightAreaEnd, sequenceLength);

    auto updateSignal = &QSpinBox::editingFinished;
    connect(sbLeftAreaStart, updateSignal, [this] { sbLeftAreaEnd->setMinimum(sbLeftAreaStart->value()); });
    connect(sbLeftAreaEnd, updateSignal, [this] { sbLeftAreaStart->setMaximum(sbLeftAreaEnd->value()); });
    connect(sbRightAreaStart, updateSignal, [this] { sbRightAreaEnd->setMinimum(sbRightAreaStart->value()); });
    connect(sbRightAreaEnd, updateSignal, [this] { sbRightAreaStart->setMaximum(sbRightAreaEnd->value()); });
}

QString PCRPrimerDesignForDNAAssemblyOPWidget::getSelectedSequence() const {
    auto selecteItems = twGeneratedSequences->selectedItems();
    CHECK(!selecteItems.isEmpty(), QString());
    SAFE_POINT(selecteItems.size() == 1, L10N::nullPointerError("Unexpected selection"), QString());

    auto sequence = selecteItems.first()->text();
    SAFE_POINT(sequence.size() == 8, L10N::nullPointerError("Unexpected sequence size"), QString());

    return sequence;
}
}
